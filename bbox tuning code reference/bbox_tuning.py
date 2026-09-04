import os
import cv2
import math
import random
import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
from scipy import stats
from collections import OrderedDict

from sklearn.metrics.pairwise import cosine_similarity

import warnings

import matplotlib as mpl
from matplotlib import pyplot as plt
import plotly.express as px

from skimage.metrics import peak_signal_noise_ratio as psnr
#from skimage.feature import peak_local_max
#from skimage.metrics import structural_similarity as ssim
#from skimage.metrics import hausdorff_distance

from skimage.filters import unsharp_mask

import trackpy as tp

import multiprocessing
from joblib import Parallel, delayed

from utils import *

plt.rcParams.update({'font.size': 18})

import logging
FS_logger = logging.getLogger('find_and_supress') #FS_logger = logging.getLogger(__name__)
IFS_logger = logging.getLogger('iterative_find_and_supress')
video_logger = logging.getLogger('annotate_objects_in_video')

def warm_start_parallel_workers(search_config, n_jobs=None):
    if not search_config.get('parallel_processing', False):
        return  # Do nothing if parallelization is off
    
    try:
        #print("Warming up parallel workers...")
        # Launch lightweight dummy tasks to warm up joblib workers
        if n_jobs is None:
            n_jobs = max(1, multiprocessing.cpu_count() // 2)
            
        Parallel(n_jobs=n_jobs)(
            delayed(lambda x: x**2)(i) for i in range(n_jobs)
        )
        #print("Parallel workers are ready.")
    except Exception as e:
        print(f"Warning: Failed to warm-start workers: {e}")

def run_balanced_parallel(jobs, max_jobs=None):
    """
    Execute jobs in parallel using joblib with optional job count control.

    Parameters
    ----------
    jobs : list
        List of delayed joblib jobs to execute.
    max_jobs : int or None
        Number of parallel processes to use. If None, uses half the CPU cores.

    Returns
    -------
    list
        List of job results.
    """
    if max_jobs is None:
        max_jobs = max(1, multiprocessing.cpu_count() // 2)

    return Parallel(n_jobs=max_jobs, batch_size=1, backend='loky')(jobs)
    

def thin_plate_spline(arr: np.ndarray, scale_by: int = 10, show_plot: bool = False) -> np.ndarray:
    """
    Apply thin plate spline interpolation to a 1D array and optionally visualize the result.

    Parameters:
    -----------
    arr : np.ndarray
        1D input array to be interpolated.
    scale_by : int, optional (default=10)
        Factor by which to increase the resolution of the output.
    show_plot : bool, optional (default=False)
        If True, displays a plot comparing the original and interpolated data.

    Returns:
    --------
    np.ndarray
        Interpolated array with higher resolution.
    """
    if arr is None or arr.ndim != 1:
        raise ValueError("Input must be a non-empty 1D NumPy array.")

    x_original = np.arange(len(arr))
    rbf_func = Rbf(x_original, arr, function='thin_plate')

    # Generate interpolated x values
    x_interp = np.linspace(0, len(arr) - 1, ((len(arr) - 1) * scale_by) + 1)
    arr_interp = rbf_func(x_interp)

    if show_plot:
        plt.figure(figsize=(10, 4))
        plt.plot(x_original, arr, 'o-', label='Original', color='blue')
        plt.plot(x_interp, arr_interp, '-', label='Thin Plate Spline', color='red')
        plt.legend()
        plt.title('Thin Plate Spline Interpolation')
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.grid(True)
        plt.show()

    return arr_interp


def create_masked_slice(mean_slice):
    mask = np.ones(len(mean_slice))
    mask[:3] = 0
    mask[-3:] = 0
    mean_slice = mean_slice*mask
    
    return mean_slice


def finding_edges(
    im_crop,
    edge_relThresh=0.15,
    tps_scale_by=10,
    center_around='Center',
    sharpen_image=True
):
    """
    Estimate object boundaries from a cropped image using thin-plate spline smoothing.

    Parameters:
    -----------
    im_crop : np.ndarray
        Cropped grayscale image.
    edge_relThresh : float
        Relative threshold to determine edge based on intensity.
    tps_scale_by : int
        Upsampling factor for smooth spline curves.
    center_around : str
        Strategy for vertical/horizontal line selection ('Peak' or 'Center').
    sharpen_image : bool
        Whether to apply unsharp masking before analysis.

    Returns:
    --------
    x_tl_new, x_br_new, y_tl_new, y_br_new : float
        Refined coordinates of the bounding box edges.
    """
    if im_crop.shape[0] <= 5 or im_crop.shape[1] <= 5:
        #warnings.warn("Crop too small in Finding edges", RuntimeWarning)
        return 0, 1, 0, 1

    # Optional sharpening (improves boundary clarity in soft-focus images)
    if sharpen_image:
        im_crop = unsharp_mask(im_crop, radius=max(im_crop.shape), amount=25)
        im_crop = im_crop - np.nanmin(im_crop)
        im_crop = im_crop * (255.0 / im_crop.max())

    # Estimate center of analysis based on Peak or Center
    if center_around == 'Peak':
        y_slice_masked = create_masked_slice(np.mean(im_crop, axis=1))
        x_slice_masked = create_masked_slice(np.mean(im_crop, axis=0))
        y_slice_center = np.argmax(y_slice_masked)
        x_slice_center = np.argmax(x_slice_masked)
    elif center_around == 'Center':
        y_slice_center = im_crop.shape[0] // 2
        x_slice_center = im_crop.shape[1] // 2
    else:
        raise ValueError("center_around must be either 'Peak' or 'Center'")

    # Extract vertical and horizontal line slices through center
    y_slice = im_crop[:, max(0, x_slice_center - 3):x_slice_center + 2]
    x_slice = im_crop[max(0, y_slice_center - 3):y_slice_center + 2, :]

    if y_slice.size == 0 or y_slice.max() == 0:
        return 0, 1, 0, 1
    if x_slice.size == 0 or x_slice.max() == 0:
        return 0, 1, 0, 1

    y_slice = np.mean(y_slice, axis=1)
    x_slice = np.mean(x_slice, axis=0)
    
    # If either slice contains only NaNs or has no valid range, return early
    if np.isnan(y_slice).all() or np.isnan(x_slice).all():
        return 0, 1, 0, 1

    # Compute relative threshold
    image_threshold = ((im_crop.max() - im_crop.min()) * edge_relThresh) + im_crop.min()

    # Apply thin-plate spline smoothing
    x_spline = thin_plate_spline(arr=x_slice, scale_by=tps_scale_by)
    y_spline = thin_plate_spline(arr=y_slice, scale_by=tps_scale_by)

    # If edges are detectable, return cropped bounds
    if image_threshold < x_spline.max() and image_threshold < y_spline.max():
        # Find the first and last x index where the spline exceeds the threshold 
        x_tl_new = np.argmax(x_spline > image_threshold) / tps_scale_by
        x_br_new = (len(x_spline) - np.argmax(x_spline[::-1] > image_threshold) - 1) / tps_scale_by

        # Find the first and last y index where the spline exceeds the threshold 
        y_tl_new = np.argmax(y_spline > image_threshold) / tps_scale_by
        y_br_new = (len(y_spline) - np.argmax(y_spline[::-1] > image_threshold) - 1) / tps_scale_by

        return x_tl_new, x_br_new, y_tl_new, y_br_new
    else:
        return 0, 1, 0, 1


#It uses match_template_in_image so prolly has to stay here for now
def search_location_with_correlation(
    frame_to_use: np.ndarray,
    previousFrame: np.ndarray,
    df_row: pd.Series,
    search_range: int = 10,
    minTemplateSize: int = 7
) -> tuple:
    """
    Use template matching to locate an object in the current frame based on its previous location.

    Parameters:
    -----------
    frame_to_use : np.ndarray
        The current video frame (grayscale).
    previousFrame : np.ndarray
        The previous video frame (grayscale).
    df_row : pd.Series
        DataFrame row with bounding box details: ['x top left', 'y top left', 'width', 'height'].
    search_range : int
        Number of pixels around the object to search in the current frame.
    minTemplateSize : int
        Minimum template dimension to ensure robust matching.

    Returns:
    --------
    tuple : (y_tl, x_tl)
        Updated top-left coordinates of the object in the current frame.
    """
    # Determine extra padding if the template is too small
    template_offset = 0
    if df_row['height'] < minTemplateSize:
        template_offset = (minTemplateSize - df_row['height']) // 2
    elif df_row['width'] < minTemplateSize:
        template_offset = (minTemplateSize - df_row['width']) // 2

    # Extract the object template from the previous frame
    object_template = row_crop(previousFrame, df_row, offset=template_offset)

    # Pad the current frame to allow searching beyond the original bounds
    image_padded = np.pad(frame_to_use, search_range, mode='constant')

    # Adjust coordinates in the DataFrame row to account for padding
    df_row_padded = df_row.copy()
    df_row_padded[['y top left', 'x top left']] += search_range

    # Crop the padded image around the original object's vicinity
    search_region = row_crop(image_padded, df_row_padded, offset=search_range)

    # Match the template in the cropped region
    search_config = None
    df_correlation_peak = match_template_in_image(
        search_region, object_template, float('nan'), df_row, search_config, detection_type="single"
    )

    # Adjust the match result to account for template offset
    df_correlation_peak['y top left'] += template_offset
    df_correlation_peak['x top left'] += template_offset

    # Update the original coordinates by shifting them based on match location
    df_row['y top left'] += df_correlation_peak['y top left'].item() - search_range
    df_row['x top left'] += df_correlation_peak['x top left'].item() - search_range

    return tuple(df_row[['y top left', 'x top left']].tolist())
    
   
def dilated_local_max(res, threshold, window_size=5):
    """
    Find local maxima in a 2D correlation map using dilation and thresholding.

    Parameters
    ----------
    res : np.ndarray
        A 2D array (usually the result of cv2.matchTemplate) representing match scores.
    threshold : float
        Minimum score required for a point to be considered a local maximum.
    window_size : int, optional
        Size of the window used for morphological dilation. Determines the minimum distance
        between detected peaks. Must be an odd positive integer. Default is 3.

    Returns
    -------
    np.ndarray
        An array of shape (N, 2), where each row is the (y, x) coordinate of a local maximum
        that meets or exceeds the threshold.
    
    Notes
    -----
    This function performs non-maximum suppression by comparing each value in `res`
    to its local neighborhood defined by `window_size`. It is significantly faster than
    `skimage.feature.peak_local_max` and does not sort or rank by score.
    """
    dilated = cv2.dilate(res, np.ones((window_size, window_size), np.uint8))
    mask = (res == dilated) & (res >= threshold)
    return np.argwhere(mask)


def shortlist_augmentations(template, flips, rotation_map, similarity_threshold=0.95):
    augmented_templates = []
    augmented_params = []

    for flip in flips:
        flipped = cv2.flip(template, flip) if flip is not None else template
        for rot_code in rotation_map:
            rotated = rotation_map[rot_code](flipped)
            resized = cv2.resize(rotated, template.shape[::-1], interpolation=cv2.INTER_AREA)
            resized = flipped.astype(np.float32).ravel()

            augmented_templates.append(resized)
            augmented_params.append((flip, rot_code))

    augmented_templates = np.array(augmented_templates)
    similarity_matrix = cosine_similarity(augmented_templates)

    # Greedy selection: pick one, then remove highly similar others
    selected_indices = []
    indices_to_remove = set()

    for i in range(len(augmented_templates)):
        if i in indices_to_remove:
            continue
        selected_indices.append(i)
        # Mark all similar variants as indices_to_remove
        similar = np.where(similarity_matrix[i] > similarity_threshold)[0]
        indices_to_remove.update(similar)

    # Get final selected (flip, rot_code) pairs
    selected_params = [augmented_params[i] for i in selected_indices]
    return selected_params


def match_one_combination(img, template, new_shape, flip, rot_code, row, threshold, rotation_map, min_size=4):
    if new_shape[0] < min_size or new_shape[1] < min_size:
        return []
    
    base_template = cv2.resize(template, new_shape, interpolation=cv2.INTER_AREA)
    flipped = cv2.flip(base_template, flip) if flip is not None else base_template
    rotated = rotation_map[rot_code](flipped)

    th, tw = rotated.shape
    res = cv2.matchTemplate(img, rotated, cv2.TM_CCOEFF_NORMED)
    #loc = np.where(res >= threshold)

    coordinates = dilated_local_max(res, threshold, window_size=3)

    return [{
        "x top left": int(x), "y top left": int(y),
        "width": tw, "height": th,
        "Aspect ratio": round(tw / th, 2) if th != 0 else 0,
        "frame": row["frame"], "madeBy": "Matching",
        "humanMade": False, "nms applied": False,  # more accurately yes now
        "new object": True,
        "match score": float(res[y, x])
    } for y, x in coordinates]


def match_template_in_image(img, template, threshold, row, search_config, detection_type="multiple", max_templates=1000):
    """
    Match a template within an image and return bounding box results.

    Parameters:
    -----------
    img : np.ndarray
        Image in which to search for the template.
    template : np.ndarray
        Template to match.
    threshold : float
        Correlation threshold for matching.
    row : pd.Series
        Reference row containing metadata like frame number.
    detection_type : str
        Type of detection: "single", "multiple", or "multiscale".

    Returns:
    --------
    pd.DataFrame
        Detected bounding boxes and metadata.
    """
    h, w = template.shape
    results = []

    if detection_type == "multiple":
        res = cv2.matchTemplate(img, template, cv2.TM_CCOEFF_NORMED)
        loc = np.where(res >= threshold)
        results = [
            {
                "x top left": pt[0], "y top left": pt[1],
                "width": w, "height": h,
                "Aspect ratio": round(w / h, 2) if h != 0 else 0,
                "frame": row["frame"], "madeBy": "Matching",
                "humanMade": False, "nms applied": False, "new object": True
            }
            for pt in zip(*loc[::-1])
        ]
    
    elif detection_type == "multiscale":
        flips = [None, 0]
        
        rotations = [None, 0, 1, 2]  # 0=90, 1=180, 2=270 degrees clockwise
        rotation_map = {
            None: lambda img: img,
            0: lambda img: cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE),
            1: lambda img: cv2.rotate(img, cv2.ROTATE_180),
            2: lambda img: cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)
        }
        
        # Remove highly similar augmentations
        selected_params = shortlist_augmentations(template, flips, rotation_map)
        
        scales = np.arange(0.5, 2.01, 0.25)
        new_shapes = sorted({tuple([int(scale * dim) for dim in template.shape])[::-1] for scale in scales})
        
        # Create list of tasks instead of nested for loops
        param_grid = [
            (new_shape, flip, rot)
            for new_shape in new_shapes
            for (flip, rot) in selected_params
        ]
        
        # Match combinations
        if search_config['parallel_processing']:
            jobs = [
                delayed(match_one_combination)(img, template, shape, flip, rot, row, threshold, rotation_map)
                for shape, flip, rot in param_grid
            ]
            all_results = run_balanced_parallel(jobs)
        else:
            all_results = [
                match_one_combination(img, template, shape, flip, rot, row, threshold, rotation_map)
                for shape, flip, rot in param_grid
            ]
        
        results = [item for sublist in all_results for item in sublist]
        
        # If I have too many templates, randomly select max_templates
        if len(results) > max_templates:
            results = random.sample(results, max_templates)
        #print(f"{len(results) = }")
         
    elif detection_type == "single":
        res = cv2.matchTemplate(img, template, cv2.TM_CCOEFF_NORMED)
        _, max_val, _, max_loc = cv2.minMaxLoc(res)
        results.append({
            "x top left": max_loc[0], "y top left": max_loc[1],
            "width": w, "height": h,
            "Aspect ratio": round(w / h, 2) if h != 0 else 0,
            "frame": row["frame"], "madeBy": "Matching",
            "humanMade": False, "nms applied": False, "new object": True,
            "match score": float(max_val)
        })

    if not results:
        return pd.DataFrame(columns=[
            "x top left", "y top left", "width", "height", "Aspect ratio",
            "frame", "madeBy", "humanMade", "nms applied", "new object", "match score"
        ])
    else:
        return pd.DataFrame(results)[[
            "x top left", "y top left", "width", "height", "Aspect ratio",
            "frame", "madeBy", "humanMade", "nms applied", "new object", "match score"
        ]]


#Implementing non-maximum supression
#Make sure to test any changes made to this robustly
def nms_with_area(df, overlapThresh, areaDiff):
    """
    Perform Non-Maximum Suppression (NMS) with area-based filtering.

    This function removes overlapping bounding boxes (bboxes) based on IoU and area differences.
    Bboxes with drastically different areas are not suppressed, even if they overlap.
    
    This happens when the objects with different areas overlap
    
      Parameters:
    ----------
    df : pandas.DataFrame
        DataFrame with bbox coordinates and dimensions. Columns required:
        'x top left', 'y top left', 'width', 'height'.
    overlapThresh : float
        IoU threshold (0 to 1) for suppressing overlapping bboxes.
    areaDiff : float
        Maximum allowable area (0 to 1) ratio between bboxes for suppression.

    Returns:
    -------
    pandas.DataFrame
        Filtered DataFrame with non-suppressed bboxes.
    """
    
    # if there are no boxes, return an empty list
    if len(df) == 0:
        return []
    
    #Drop any duplicates
    df = df.drop_duplicates().reset_index(drop=True)
    
    #df = df.sort_values(by=['madeBy'], ascending=[True]).reset_index(drop=True)
    
    # initialize the list of picked indexes	
    pick = df[df['nms applied']].index.to_list() #[]
    #if len(pick) > 0:
    #    print(f'Total Items = {len(df)}, Already Found: {len(pick)}')
    #print('Already Found:', len(pick), pick)
    #temp = df[df['nms applied']].index.to_list() #
    #print(f"{len(temp)}, {temp} = ")
    
    # grab the coordinates of the bounding boxes
    x1 = df['x top left'].astype("float").to_list()
    y1 = df['y top left'].astype("float").to_list()
    x2 = (df['x top left']+df['width']).astype("float").to_list() 
    y2 = (df['y top left']+df['height']).astype("float").to_list() 
    #a
    
    width = df['width'].astype("float").to_list()
    height = df['height'].astype("float").to_list()
    #a
    
    # compute the area of the bounding boxes
    area = [width[i]*height[i] for i in range(len(width))] #(x2 - x1 + 1) * (y2 - y1 + 1)
    #print(area)
    
    # sort the bounding boxes by the bottom-right y-coordinate of the bounding box
    #The df is already sorted with top-left x-coordinate. Keeping it that way
    #idxs are the true indices of the df (Note that these aren't always sorted)
    #idxs = np.argsort(y2)
    idxs = df.index.to_list() #np.argsort(x1)
    #print(f"{len(idxs) = }, Unique items = {len(list(set(idxs)))}")
    #temp = 0
    #a
    
    # keep looping while some indexes still remain in the indexes list
    #Only keeping the last box in a vicinity
    #If an overlap is found, suppress that location of index and don't loop over it
    #Hence, some indices locations are actively removed from while
    #This contributes to actively speeding up the loops
    #Note 1: Pick keeps the indices to keep in the filtered df after the for loop
    #Note 2: Suppress removes some index locations to actively supress the high overlap locations of the indices
    #save_flag = True
    while len(idxs) > 0:
        # grab the last index in the indexes list and add the
        # index value to the list of picked indexes

        #Last is not the true index of the df
        #It is just a number representing where this element is in the list
        #Don't use it as an index of df
        last = len(idxs) - 1
        #i is the right index that can be used to select a row from the df
        i = idxs[last]

        #Save the index if it is not in the pick list already
        if i not in pick:
            pick.append(i)
        suppress = [last]

        #temp += 1

        #print(f"Keeping: {i}", df.loc[i].to_dict())
        #print(df.loc[i].to_dict())
        #print(nms_applied[last])

        #if nms_applied[last]:
        #    continue
        
        # find the largest (x, y) coordinates for the start of the bounding box and the 
        # smallest (x, y) coordinates for the end of the bounding box loop over all indexes 
        # in the indexes list
        #As len(pos) = last, the max value of pos is 1 less than last
        #Hence, i != j
        for pos in range(0, last):
            #j is one of index of the df. It can be used to select a row
            j = idxs[pos]
            
            #print(pos, j)
            #a
            
            #If areas are drastically different, don't suppress the index
            #This is probably a big object overlapping a smaller one
            if area[i] < areaDiff*area[j]:
                continue
            if area[j] < areaDiff*area[i]:
                continue
                
            # find the largest (x, y) coordinates for the start of the bounding box and 
            # the smallest (x, y) coordinates for the end of the bounding box
            xx1 = max(x1[i], x1[j])
            yy1 = max(y1[i], y1[j])
            xx2 = min(x2[i], x2[j])
            yy2 = min(y2[i], y2[j])
            
            # compute the width and height of the bounding box
            w = max(0, xx2 - xx1 + 1)
            h = max(0, yy2 - yy1 + 1)
            # compute the ratio of overlap between the computed
            # bounding box and the bounding box in the area list
            overlap = float(w * h) / area[j]
            
            # if there is sufficient overlap, suppress the location of the index
            #Note this process is not removing the row of that bbox
            #It is only removing the location of that index so that it does not end up in the pick list
            if overlap > overlapThresh:
                suppress.append(pos)
                #print(f'Supressing {pos} based on {last}')
            #else: #A lot of these are bad aspect ratios. 
            #We are keeping it for now as it might still be useful for object detection
                #if abs(x1[i] - x1[j]) == 1 and abs(y1[i] - y1[j]) == 1:
                #    print(i, "i: ", df.loc[i].to_dict())
                #    print(j, "j: ", df.loc[j].to_dict())
                #    print(w, h, area[j], overlap)
                #    if save_flag:
                #        df.to_evxcel(f'temp_{len(df)}.xlsx', index=False)
                #        save_flag = False
                        
        
        # delete all indexes from the index list that are in the suppression list
        idxs = np.delete(idxs, suppress)
    
    #print(f"{temp = }")
    #print(pick)
    #print('Matching:', len(df), df.index.to_list()[:50])
    #print('Picked: ', len(pick), sorted(pick))
    df = df[df.index.isin(pick)]
    #print('NMS:', len(df), df.index.to_list(), '\n')
    ###print('x top left:', df['x top left'].to_list(), '\n')
    
    # Mark remaining new rows as processed
    #df['nms applied'] = True
    
    return df
    

def detect_and_refine_boundary(image, x_tl, y_tl, w, h, pixel_pad=5, edge_relThresh=0.15, 
                                tps_scale_by=10, temporal_detection=False, visualize=False): #x_topLeft #width_multiplier, 
                                #ap_relThresh=0.75, 
    """
    This code
     - uses initial co-ordinates as inputs
     - Obtain more locations by shaking the box around the initial location
     - Make crops at those locations and find edges in those locations/crops
     - Returns the location for the found object
     - Provide geometric and intensity center for the found object
    """
    
    #Return No detection if the inputs are nan
    if math.isnan(x_tl) or math.isnan(y_tl) or math.isnan(w) or math.isnan(h):
        return 0, 0, 0, 0, 'No detection', False, False
    
    #Extracting single channel (green)
    if len(image.shape) == 3:
        image = image[:,:,1] #Works for Ardekani, Hydrogel, and KIT data
        image = pre_processing(image)
    
    #Adding a zero-pad around the image
    image_padded = np.pad(image, pixel_pad, 'constant')
    #print(image.shape, image_padded.shape)

    hCropStart = max(1, math.floor(x_tl)) #math.floor(row['x_im']) # - pixel_pad)
    hCropEnd = math.ceil(x_tl + w) + 2*pixel_pad
    vCropStart = max(1, math.floor(y_tl)) # - pixel_pad)
    vCropEnd = math.ceil(y_tl + h) + 2*pixel_pad
    #print(y_tl, h, x_tl, w)
    #print(vCropStart, vCropEnd, hCropStart, hCropEnd)
    
    im_crop = padded_crop(image_padded, vCropStart, vCropEnd, hCropStart, hCropEnd) #, ap_relThresh
    #display_image(im_crop, style = "plotly") #
    
    edge_Thresh = edge_relThresh
    x_tl_new, x_br_new, y_tl_new, y_br_new = finding_edges(im_crop, edge_Thresh, tps_scale_by, 
                                                            sharpen_image=False)
    #print(x_tl_new, x_br_new, y_tl_new, y_br_new)
    
    #Redetect if there is a finding_edges failure or object width or height are too small
    if (x_tl_new == 0 and x_br_new == 1) or (x_tl_new - x_br_new < 2) or (y_tl_new - y_br_new < 2):
        #print("Redection with sharpen")
        x_tl_new, x_br_new, y_tl_new, y_br_new = finding_edges(im_crop, edge_Thresh, tps_scale_by, 
                                                            sharpen_image=True)
    
    #As the objects don't disappear temporally (except on image boundaries)
    #There is a high likelihood that the object exists and the edge detection algorithm failed
    #Trying to detect the object again with a lower threshold.
    if temporal_detection:
        if x_tl_new == 0 and x_br_new == 1:
            edge_Thresh = edge_relThresh/2
            x_tl_new, x_br_new, y_tl_new, y_br_new = finding_edges(im_crop, edge_Thresh, tps_scale_by, 
                                                                sharpen_image=False)
                                                            
        
    #Adding image dimensions and removing padding
    #Padding is already a part of the big crop
    y_tl_new = max(0, y_tl_new + y_tl - pixel_pad)
    y_br_new = min(y_br_new + y_tl - pixel_pad, image.shape[0]-1)
    x_tl_new = max(0, x_tl_new + x_tl - pixel_pad)
    x_br_new = min(x_br_new + x_tl - pixel_pad, image.shape[1]-1)
    #print(y_tl_new, y_br_new, x_tl_new, x_br_new, image.shape)
        
    if y_tl_new > image.shape[0]-2 or x_tl_new > image.shape[1]-2:
        return 0, 0, 0, 0, 'No detection', False, False
        
    #print(y_tl_new, y_br_new, x_tl_new, x_br_new, image.shape)
    h_new = round(y_br_new-y_tl_new, 2)#()/tps_scale_by
    w_new = round(x_br_new-x_tl_new, 2)#()/tps_scale_by

    madeBy = 'Detected'

    #Detection can cause nearby bboxes to converge at a single location
    #These should be eliminated in the next NMS round
    nms_applied = False
    
    return y_tl_new, x_tl_new, h_new, w_new, madeBy, edge_Thresh, nms_applied #x_tl, y_tl, w, h


def refine_boundary(currentFrame, row, edge_detection_config):
    matched_object = row_crop(currentFrame, row, offset=0)
    y_tl, x_tl, h, w = row[['y top left', 'x top left', 'height', 'width']]

    y, x, h, w, madeBy, edgeThreshold, nmsApplied = detect_and_refine_boundary(
        currentFrame, x_tl, y_tl, w, h, pixel_pad=edge_detection_config['pixel_pad'],
        edge_relThresh=edge_detection_config['edge_relThresh']
    )

    row_data = row.to_dict()
    row_data.update({ 'y top left': y, 'x top left': x, 'height': h, 'width': w, 
                     'madeBy': madeBy, 'edgeThreshold': edgeThreshold, 
                     'nms applied': nmsApplied
                    })

    if madeBy != 'No detection':
        detected_object = cv2.resize(row_crop(currentFrame, row_data, offset=0),
            dsize=matched_object.shape[::-1], interpolation=cv2.INTER_CUBIC
                                    )
        
        row_data['Detection Similarity'] = cv2.matchTemplate(matched_object, detected_object, 
                                                             cv2.TM_CCOEFF_NORMED)[0][0]
        
        row_data['PSNR'] = psnr(matched_object, detected_object)

    return row_data
    

def find_and_supress(df, currentFrame, match_threshold, max_templates, edge_detection_config, search_config, 
                    postprocess_config):
    df_new_object = df[df['new object']].copy()
    if len(df_new_object) == 0:
        return df
    
    row2 = df_new_object.sample(frac=1).iloc[0]
    currentTemplate = row_crop(currentFrame, row2, offset=0)
    
    start_time = time.time()

    #Mask the image based on existing annotations
    mask = create_mask_from_annotations(df, currentFrame.shape, pad=0)
    maskedFrame = apply_mask_to_image(currentFrame, mask, mode='black') #change to none for no mask
    
    #maskedFrame
    
    df_found = match_template_in_image(maskedFrame, currentTemplate, match_threshold, row2, search_config, 
                                       detection_type="multiscale", max_templates=max_templates)
    FS_logger.info(f"Match Template|||{time.time()-start_time:.2f}|||{len(df_found)}")
    
    #Concat with previous df, remove duplicates, and keep only the new objects
    if len(df_found) > 0:
        #Make a temporary df with nms applied
        df_temp = df.copy() 
        df_temp['nms applied'] = True
        
        df_found = pd.concat([df_temp, df_found], ignore_index=True)
        df_found = df_found.sort_values(by=['humanMade', 'nms applied', 'x top left', 'y top left'],
                                          ascending=[True, True, True, True]).reset_index(drop=True)

        start_time = time.time()
        df_found = nms_with_area(
            df_found, postprocess_config['NMSoverlapThresh'], postprocess_config['NMSareaDiff']
        )
        
        ############ I am not effectively removing the df_temp elements... I should remove it and merge it larer
        df_found = df_found[df_found['new object']]
        FS_logger.info(f"NMS after Matching|||{time.time()-start_time:.2f}|||{len(df_found)}")
        
        if len(df_found) > 0:
            df_found = df_found[[
                'x top left', 'y top left', 'width', 'height', 'frame',
                'madeBy', 'humanMade', 'nms applied', 'new object'
            ]]
    
    #df_found.to_excel(os.path.join('output_files', 'IMG_1301_frames30583-30883_firstFrame_b.xlsx'), index=False)
    
    #Refine boundaries
    start_time = time.time()
    if search_config['parallel_processing'] and len(df_found) > 1:
        jobs = [delayed(refine_boundary)(currentFrame, row, edge_detection_config)
            for _, row in df_found.iterrows()
               ]
        refined_rows = run_balanced_parallel(jobs)
        df_found = pd.DataFrame(refined_rows)
        
    else:
        refined_rows = [
            refine_boundary(currentFrame, row, edge_detection_config)
            for _, row in df_found.iterrows()
        ]
        df_found = pd.DataFrame(refined_rows)
    FS_logger.info(f"Refine boundaries (No Drop)|||{time.time()-start_time:.2f}|||{len(df_found)}")
    #df_found.to_excel(os.path.join('output_files', 'IMG_1301_frames30583-30883_firstFrame_c.xlsx'), index=False)
        
    if len(df_found) == 0:
        return df
    
    df_found = df_found[df_found['width'] > 2]
    df_found = df_found[df_found['height'] > 2]
    
    df_found['Aspect ratio'] = [round(item, 2) for item in df_found['width'] / df_found['height']]
    
    if len(df_found) > 0:
        if 'PSNR' in df_found.columns:
            df_found = df_found[df_found['PSNR'] < postprocess_config['PSNR']]

    ar_min, ar_max = postprocess_config['ar_range']
    if len(df_found) > 0:
        df_found = df_found[(df_found['Aspect ratio'] > ar_min) & (df_found['Aspect ratio'] < ar_max)]

    FS_logger.info(f"Edge Detection and Applying Filters|||{time.time()-start_time:.2f}|||{len(df_found)}")
    
    #if new found objects remaining
    if len(df_found) > 0:
        start_time = time.time()
        df_found['objectID'] = df_found.index
        df_found = nms_with_area(df_found, postprocess_config['NMSoverlapThresh'], postprocess_config['NMSareaDiff']
                    ).reset_index(drop=True)
        FS_logger.info(f"NMS after Filters and Edge Detection|||{time.time()-start_time:.2f}|||{len(df_found)}")
        
        #Removing unstable annotations from the newly found ones
        #Keeping this is useful as it removes some of the outliers even in the first frame
        #This function is especially important when we are not doing any peak-based detection
        start_time = time.time()
        df_found = remove_unstable(currentFrame, df_found, edge_detection_config, parallel_processing=search_config['parallel_processing'])
        FS_logger.info(f"Unstable removed|||{time.time()-start_time:.2f}|||{len(df_found)}")
    
    #df_found.to_excel(os.path.join('output_files', 'IMG_1301_frames30583-30883_firstFrame_d.xlsx'), index=False)
    
    #if new found objects still remaining
    if len(df_found) > 0:
        df = pd.concat([df, df_found], ignore_index=True)
        FS_logger.info(f"Concat. Before final NMS|||{time.time()-start_time:.2f}|||{len(df)}")
        df = df.sort_values(by=['humanMade', 'nms applied', 'x top left', 'y top left'],
                            ascending=[True, True, True, True]).reset_index(drop=True)
    
        df = nms_with_area(df, postprocess_config['NMSoverlapThresh'], postprocess_config['NMSareaDiff']).reset_index(drop=True)
        FS_logger.info(f"After final NMS|||{time.time()-start_time:.2f}|||{len(df)}")
        
        df['nms applied'] = True
        df['objectID'] = df.index
        #FS_logger.info(f"Final Merge|||{time.time()-start_time:.2f}|||{len(df)}")
    
    return df
    
    
def refine_human_annotation(
    df_new_annotation: pd.DataFrame,
    df: pd.DataFrame,
    currentFrame: np.ndarray,
    edge_detection_config: dict,
    search_config: dict,
    tabular_save_path: str = None,
) -> pd.DataFrame:
    """
    Iteratively detect similar objects starting from an initial human annotation.

    Parameters:
    -----------
    df_new_annotation : pd.DataFrame
        New annotation(s) (typically from a human).
    currentFrame : np.ndarray
        Image frame used for detection.
    edge_detection_config : dict
        Parameters for detection (pixel_pad, edge threshold, etc).
    search_config : dict
        Parameters for running the loops (max_iterations, iteration_memory, etc.)
    tabular_save_path : str or None
        Optional path to save result after processing.

    Returns:
    --------
    pd.DataFrame
        Final DataFrame with refined annotations.
    """
    if len(df_new_annotation) != 1:
        print("More than 1 annotation provided. Please check the workflow.")
        return df
    
    # Extract initial box
    y_tl, x_tl, h, w = df_new_annotation.loc[0, ['y top left', 'x top left', 'height', 'width']]
    base_thresh = edge_detection_config['edge_relThresh']
    
    # List of edge thresholds to try, capped at 1.0
    thresh_scales = [1, 2, 1/2] #3, 1/3 #[1] #
    #stable_detection_found = False
    df_return_annotation = None
    
    for scale in thresh_scales:
        edge_thresh = min(scale * base_thresh, 1.0) #Max threshold should be 1
        
        #if scale != 1:
        #    print(f"Changing edge threshold by {scale}x")

        # Run edge detection
        y, x, h_new, w_new, madeBy, edgeThreshold, nms_applied = detect_and_refine_boundary(
            currentFrame, x_tl, y_tl, w, h, pixel_pad=0,
            edge_relThresh=edge_thresh
        )

        # Abort early if objects are too small
        if h_new < 2 or w_new < 2:
            print(f"Low dimension detected {h_new = } and {w_new = }")
            continue

        df_test_annotation = df_new_annotation.copy()
        df_test_annotation.loc[0, ['y top left', 'x top left', 'height', 'width',
                                   'madeBy', 'edgeThreshold', 'nms applied']] = \
            y, x, h_new, w_new, madeBy, edgeThreshold, nms_applied

        #Adding objectID
        df_test_annotation['objectID'] = df_test_annotation.index

        #Save the default annotation if it is None. It comes from the first successful detection
        #This will be used if all the detections turn out to be unstable
        #if scale == 1:
        if df_return_annotation is None:
            df_return_annotation = df_test_annotation.copy()

        # Test for stability
        df_stable = remove_unstable(currentFrame, df_test_annotation, edge_detection_config, parallel_processing=search_config['parallel_processing'])
        
        #Update the default annotation if the annotation is stable
        if len(df_stable) > 0:
            #stable_detection_found = True
            df_return_annotation = df_test_annotation.copy()
            break
    
    #If the code breaks early, return the original df
    if df_return_annotation is None:
        print("The object is too small to refine boundaries")
        return df
        
    # Finalize and return
    df_return_annotation = df_return_annotation.round({'x top left': 2, 'y top left': 2, 
                                                 'height': 2, 'width': 2})

    #We trust human visual system more than algorithms. So, the bbox they made should stay after refinement.
    df_return_annotation['nms applied'] = True
    
    df = pd.concat([df_return_annotation, df], ignore_index=True)
    df = df.sort_values(by=['x top left', 'y top left']).reset_index(drop=True)
    df['objectID'] = df.index

    # Save file in the specified format
    if tabular_save_path:
        _, extension = os.path.splitext(tabular_save_path)
        extension = extension[1:]
        if extension.lower() == 'parquet':
            df.to_parquet(tabular_save_path, index=False, engine='pyarrow')
        elif extension.lower() == 'xlsx':
            df.to_excel(tabular_save_path, index=False)
        else:
            raise ValueError(f"Unsupported extension '{extension}'. Use 'parquet' or 'excel'.")
            

    return df
    
    
def iterative_find_and_supress(
    df: pd.DataFrame,
    currentFrame: np.ndarray,
    matching_config: dict,
    edge_detection_config: dict,
    search_config: dict,
    postprocess_config: dict,
    tabular_save_path: str = None,
) -> pd.DataFrame:
    """
    Iteratively detect similar objects starting from an initial human annotation.

    Parameters:
    -----------
    df : pd.DataFrame
        Existing DataFrame of all found annotations.
    currentFrame : np.ndarray
        Image frame used for detection.
    matching_config: dict
        Parameters for template matching (similarity threshold, adaptive threshold change scale)
    edge_detection_config : dict
        Parameters for detection (pixel_pad, edge threshold, etc).
    search_config : dict
        Parameters for running the loops (max_iterations, iteration_memory, etc.)
    postprocess_config : dict
        Parameters for cleaning (NMS threshold, AR range, etc).
    tabular_save_path : str or None
        Optional path to save result after processing.

    Returns:
    --------
    pd.DataFrame
        Final DataFrame with refined and matched object annotations.
    """

    print('Using new object for assisted annotations')
    
    # Track annotation count over recent iterations for early stopping
    top_n = [0] * search_config['iteration_memory']
    match_threshold = matching_config['Threshold']
    max_templates = matching_config['max_templates']
    match_threshold_changed = False
        
    df_temp = df.copy()
    
    #This for loop will be run if max_iterations is at least 1
    for i in range(1, search_config['max_iterations']+1):
        start_time = time.time()
        df_len = len(df)
        
        # Attempt to detect and refine new objects
        df = find_and_supress(
            df, currentFrame, match_threshold, max_templates, edge_detection_config, search_config, postprocess_config
        )
        
        delta_time = round(time.time() - start_time, 2)
        FS_logger.info(f"Iteration {i}|||{delta_time:.2f}|||{len(df)}")
        
        print(f"Objects found after iteration {i} = {len(df)}. Δitem = {len(df) - df_len} with Δt = {delta_time}s")
        #print(f"{len(df[df['humanMade']])}|||{i}|||{delta_time:.2f}|||{len(df) - len(df_temp)}|||{len(df)}")
        IFS_logger.info(f"{len(df[df['humanMade']])}|||{i}|||{delta_time:.2f}|||{len(df) - df_len}|||{len(df)}")
        
        #Decrease the match threshold if the object does not lead to any new annotation after the first search
        if i == 1 and len(df) <= len(df_temp):
            match_threshold = matching_config['Threshold']*matching_config['adaptive_scale']
            match_threshold_changed = True
            print(f"Decreasing the similarity threshold by {matching_config['adaptive_scale']}x... New threshold {round(match_threshold, 2)}")
            #print(f"No new objects found with this annotation... Breaking the loop")
            #break
            #print("Stopping the search")
            #break
        #If new objects are found in the first iteration, bring the similarity threshold back to the original value
        if i == 2 and len(df) <= len(df_temp):
            print(f"New objects found with reduced similarity threshold... Breaking the loop")
            break
        #if i == 2 and len(df) > len(df_temp) and match_threshold_changed:
        #    match_threshold = matching_config['Threshold']
        #    print(f"New objects found... Changing the similarity threshold back to {match_threshold}")
            
            
        # Check if new objects were added
        if len(df) > min(top_n):
            top_n.append(len(df))
            top_n = top_n[1:]  # Keep the sliding window
        else:
            print(f"No new objects found in the last {search_config['iteration_memory']} iterations")
            print("Stopping the search")
            break

    #Round the dimension columns
    df = df.round({'x top left': 2, 'y top left': 2, 'height': 2, 'width': 2})
    
    #Sort and reset objectID for downstream logic
    df = df.sort_values(by=['x top left', 'y top left']).reset_index(drop=True)
    df['objectID'] = df.index
    
    
    # Save file in the specified format
    if tabular_save_path:
        if extension.lower() == 'parquet':
            df.to_parquet(tabular_save_path, index=False, engine='pyarrow')
        elif extension.lower() == 'xlsx':
            df.to_excel(tabular_save_path, index=False)
        else:
            raise ValueError(f"Unsupported extension '{extension}'. Use 'parquet' or 'excel'.")
            
    print(f"Ready for the next annotation\n")
    
    return df
    
    
def ajdust_small_object_ar(row, ar_range, size_threshold=5, range_expansion=0.25):
    ar = row['Aspect ratio']
    h, w = row['height'], row['width']

    if h == 0 or w == 0:
        return False  # invalid size → reject

    min_ar, max_ar = ar_range
    
    if h < size_threshold or w < size_threshold:
        min_adj = min_ar * (1 - range_expansion)
        max_adj = max_ar * (1 + range_expansion)
    else:
        min_adj = min_ar
        max_adj = max_ar

    return min_adj < ar < max_adj
    

def process_object_in_video(row, currentFrame, previousFrames, currentFrameID, video_annotate_config, 
                            edge_detection_config, postprocess_config):
    min_size_change = min(video_annotate_config['size_change_limit'])
    max_size_change = max(video_annotate_config['size_change_limit'])
    ar_min, ar_max = postprocess_config['ar_range']

    h, w, last_used_edgeThreshold = row[['height', 'width', 'edgeThreshold']]
    #previousFrame = previousFrames[row['frame']]
    """
    if 3 in previousFrames.keys():
        #print(f"{previousFrames.keys()}")
        print(f"{currentFrameID = }")
        print(f"{row = }")
    """
    
    #previousFrame = previousFrames[row['last_detected_frame']] #This fails for matched
    previousFrame = previousFrames[row['frame']]
    
    if video_annotate_config['initial_search_by'] == 'matching':
        y_tl, x_tl = search_location_with_correlation(currentFrame, previousFrame, row, video_annotate_config['search_range'])
    elif video_annotate_config['initial_search_by'] == 'previous_location':
        y_tl, x_tl = row[['y top left', 'x top left']]
    else:
        video_annotate_config['search_range'] = 10
        y_tl, x_tl = search_location_with_correlation(currentFrame, previousFrame, row, video_annotate_config['search_range'])

    y, x, h, w, madeBy, edgeThreshold_to_log, nmsApplied = detect_and_refine_boundary(
        currentFrame, x_tl, y_tl, w, h,
        pixel_pad=edge_detection_config['pixel_pad'],
        edge_relThresh=last_used_edgeThreshold,
        temporal_detection=True
    )

    #Some objects are under high noise or have a low contrast. 
    #Using a pad here leads to outliers in edge detection
    #If the width or height changes drastically, don't use the pad
    #Alternatively, use match only [This might be more susceptible to drift]
    if (
        h - row['height'] >= edge_detection_config['pixel_pad'] or 
        w - row['width'] >= edge_detection_config['pixel_pad'] or
        h > max_size_change * row['height'] or
        w > max_size_change * row['width']
    ):
        y, x, h, w, madeBy, edgeThreshold_to_log, nmsApplied = detect_and_refine_boundary(
            currentFrame, x_tl, y_tl, w, h,
            pixel_pad=0,
            edge_relThresh=last_used_edgeThreshold,
            temporal_detection=True
        )

    #If objects come close, the default pixel adding can become too large
    #This will result in the new object causing a high value at the border, hence misdetecting it as an edge
    #Use the initial df row, the identified row, and the pixel pad to identify it
    #If it is identified, remove the expansion for this object for this timestep to avoid the error
    if (
        round(row['y top left'], 2) == round(y + edge_detection_config['pixel_pad'], 2) or 
        round(row['x top left'], 2) == round(x + edge_detection_config['pixel_pad'], 2) or  
        round(row['y top left'] + row['height'], 2) == round(y + h - edge_detection_config['pixel_pad'], 2) or
        round(row['x top left'] + row['width'], 2) == round(x + w - edge_detection_config['pixel_pad'], 2)
    ):
        y, x, h, w, madeBy, edgeThreshold_to_log, nmsApplied = detect_and_refine_boundary(
            currentFrame, x_tl, y_tl, w, h,
            pixel_pad=math.ceil(edge_detection_config['pixel_pad']/2),
            edge_relThresh=last_used_edgeThreshold,
            temporal_detection=True
        )

    #If the aspect ratio is out of range, try with double edge_relThresh
    #ar can be replaced by some bbox consistency criterion if I choose to include it
    #I should include it (if time allows), as large shakes in width and height are not physically correct
    aspect_ratio = round(w / h, 2) if h else 0
    if aspect_ratio < ar_min or aspect_ratio > ar_max:
        y_temp, x_temp, h_temp, w_temp, madeBy_temp, edgeThresh_temp, nmsApp_temp = detect_and_refine_boundary(
            currentFrame, x_tl, y_tl, w, h,
            pixel_pad=edge_detection_config['pixel_pad'],
            edge_relThresh=1.5 * last_used_edgeThreshold,
            temporal_detection=True
        )
        #Only update this if the resulting box isn't too small
        if h > min_size_change * row['height'] and w > min_size_change * row['width']:
            y, x, h, w, madeBy, edgeThreshold_to_log, nmsApplied = (
                y_temp, x_temp, h_temp, w_temp, madeBy_temp, edgeThresh_temp, nmsApp_temp)

    #Getting similarity between match and detect
    matched_object = row_crop(previousFrame, row, offset=0)
    detect_row = pd.Series({"x top left": x, "y top left": y, "width": w, "height": h})
    detected_object = cv2.resize(row_crop(currentFrame, detect_row, offset=0), dsize=matched_object.shape[::-1],
                                 interpolation=cv2.INTER_CUBIC)
    match_score = cv2.matchTemplate(matched_object, detected_object, cv2.TM_CCOEFF_NORMED)[0][0]

    #If detection fails, we can end up with erroneous detections (pixel threshold used 2px)
    #If so, find the use matching only
    if (h <= 2 or w <= 2 or match_score < video_annotate_config['matchDetectSimilarity'] or
        h < min_size_change * row['height'] or w < min_size_change * row['width']):
        y, x, h, w = row['y top left'], row['x top left'], row['height'], row['width']
        madeBy = 'Matching'
        edgeThreshold_to_log = last_used_edgeThreshold
        nmsApplied = False
        
    #Updaing last_detected_frame
    if madeBy == 'Matching': #Use the previous one
        last_detected_frame = row['last_detected_frame']
    elif madeBy == 'Detected': #Update this
        last_detected_frame = currentFrameID
     
    aspect_ratio = round(w / h, 2) if h else 0
    return {
        'objectID': row['objectID'], 'x top left': x, 'y top left': y, 'width': w, 'height': h,
        'Aspect ratio': aspect_ratio, 'frame': currentFrameID, 'madeBy': madeBy, #row['frame'] + 1 #currentFrameID
        'last_detected_frame': last_detected_frame,
        'humanMade': False, 'edgeThreshold': edgeThreshold_to_log, 'nms applied': nmsApplied
    }


def annotate_objects_in_video(vidcap, df, video_annotate_config, edge_detection_config, search_config, postprocess_config, #detect_object_in_video
                              save_config): #, save_frequency=1, #blur_radius, 
    
    #Updating search range if initial_search_by == 'matching' and search_range==0
    if video_annotate_config['initial_search_by'] == 'matching' and video_annotate_config['search_range']==0:
        print('initial_search_by should have a search range > 0... Using a search_range of 10')
        video_annotate_config['search_range']=10
    
    if video_annotate_config['initial_search_by'] not in ['matching', 'previous_location']:
        print('initial_search_by can only be matching/previous_location')
        print('Defaulting to initial_search_by matching with a search range of 10px')
        video_annotate_config['search_range'] = 10
    
    if video_annotate_config['last_detected_frame_memory'] < 1:
        print('last_detected_frame_memory should be >= 1')
        print('Defaulting it to 10 frames')
        video_annotate_config['last_detected_frame_memory'] = 10
    # df already contains the detection from the first frame
    df_all_annotations = df.copy()
    
    print(f"Total objects: {len(df)}")
    totalFrames = int(vidcap.get(cv2.CAP_PROP_FRAME_COUNT))
    
    if video_annotate_config['frames_to_use'] > totalFrames:
        print(f'The video has only {totalFrames} frames. Using this as the maximum number')
        video_annotate_config['frames_to_use'] = totalFrames
    
    if max(video_annotate_config['size_change_limit']) < 1:
        print('As the object boundaries can change, a max size change of > 1 is recommended')
        print('As it is < 1, changing to 10% max change per frame, leading to default value of 1.1')
        video_annotate_config['size_change_limit'][0] = 1.5
    if min(video_annotate_config['size_change_limit']) > 1:
        print('As the object can shrink, a min size change of < 1 is recommended')
        print('As it is > 1, changing to 10% max change per frame, leading to default value of 0.9')
        video_annotate_config['size_change_limit'][1] = 0.5
    #print(totalFrames)
    
    previousFrames = OrderedDict()
    previousFrames[0] = load_frame(vidcap, 0)
    
    df_all_annotations['last_detected_frame'] = 0
    last_detected_frame_memory = video_annotate_config['last_detected_frame_memory']
    
    for currentFrameID in range(1, video_annotate_config['frames_to_use']): #Don't process the 0th frame
        start_time = time.time()
        
        #Get current frame
        currentFrame = load_frame(vidcap, frame_index=currentFrameID)
        previousFrames[currentFrameID] = currentFrame
        #print(f"1: {len(previousFrames) = }")
        
        #Capping the maximum number of elements in the dict
        #previousFrames also contains the current frame as it is added
        if len(previousFrames) > last_detected_frame_memory + 1:
            previousFrames.popitem(last=False)
            
        #print(f"{previousFrames.keys() = }, {len(previousFrames) = }")
        
        #Select only the df only for the last last_detected_frame_memory frames
        df_last_detected = df_all_annotations[df_all_annotations["frame"] >= currentFrameID-last_detected_frame_memory]
        
        #Get the most recent detected object
        df_detected = df_last_detected[df_last_detected["madeBy"] == "Detected"]
        df_detected = (
            df_detected.sort_values(["objectID", "last_detected_frame"], ascending=[True, False])
            .drop_duplicates("objectID", keep="first")
            )
        detected_ids = set(df_detected["objectID"])
        
        #If detected not found, get the oldest matched object
        df_matched = df_last_detected[df_last_detected["madeBy"] == "Matching"] #Get all matched
        df_matched = df_matched[~df_matched["objectID"].isin(detected_ids)] #Dropped matched if objectID is in detected

        df_matched = (
            df_matched.sort_values(["objectID", "frame"], ascending=[True, True])
            .drop_duplicates("objectID", keep="first")
        )
        
        #Merge the df
        df_last_detected = pd.concat([df_detected, df_matched], ignore_index=True).reset_index(drop=True)
        #print(f"{len(df_detected) =}, {len(df_matched) =}, {len(df_matched)+len(df_detected) = }, {len(df_last_detected) = }")
        
        #If objectID in the merged one are not equal to the unique objectIDs in the df_last_detected, 
        #missing_ids = set(df_all_annotations['objectID'].unique()) - set(df_last_detected['objectID'].unique())
        #if missing_ids:
        #    print(f"{currentFrameID = } Missing objectIDs: {missing_ids}. Breaking the loop.")
        #    break
        
        """
        #Get the df which represents the object location in the previous frame
        df_last_detected = df_all_annotations[df_all_annotations["madeBy"] == "Detected"]
        
        #Sort by objectID and frame and drop duplicates
        df_last_detected = (df_last_detected.sort_values(["objectID", "frame"], ascending=[True, False])
                            .drop_duplicates("objectID", keep="first").reset_index(drop=True)
                        )
        """
        #Read the previous frame for template creation
        #needed_frames = df_last_detected['frame'].unique()
        #previousFrames = {fid: load_frame(vidcap, fid) for fid in needed_frames}
        
        #previousFrame = load_frame(vidcap, frame_index=currentFrameID-1)
        #print(currentFrameID, currentFrame.shape, previousFrame.shape)
        
        #Create an empty list for objects identified in the new frame
        #df_current_frame = pd.DataFrame(columns=df_last_detected.columns)
        annotations = []
        
        # loop over all templates
        #Refine boundaries
        if search_config['parallel_processing']:
            jobs = [delayed(process_object_in_video)(row, currentFrame, previousFrames, currentFrameID, 
                                     video_annotate_config, edge_detection_config, postprocess_config)
                for _, row in df_last_detected.iterrows()
                   ]
            refined_rows = run_balanced_parallel(jobs)
            df_current_frame = pd.DataFrame([r for r in refined_rows if r is not None])
            
        else:
            refined_rows = [
                process_object_in_video(row, currentFrame, previousFrames, currentFrameID, video_annotate_config,
                                edge_detection_config, postprocess_config)
                for _, row in df_last_detected.iterrows()
            ]
            df_current_frame = pd.DataFrame([r for r in refined_rows if r is not None])
        
        #print(f"{currentFrameID = } \t {len(df_current_frame) = }")
        
        #Cleaning by aspect ratio (Useful to remove false positives early for the object dataset)
        #If the object is small, the aspect ratio is a bit relaxed
        df_current_frame = df_current_frame[
            df_current_frame.apply(lambda row: ajdust_small_object_ar(row, postprocess_config['ar_range']), axis=1)
        ]
        
        #df_current_frame = df_current_frame[df_current_frame['Aspect ratio'] > 
        #                                    min(postprocess_config['ar_range'])]
        #df_current_frame = df_current_frame[df_current_frame['Aspect ratio'] < 
        #                                    max(postprocess_config['ar_range'])]
        
        df_all_annotations = pd.concat([df_all_annotations, df_current_frame]).reset_index(drop=True) #, ignore_index=True
        
        delta_time = round(time.time() - start_time, 2)
        print(f"Objects found in frame {currentFrameID}: {len(df_current_frame)}, Running Total = {len(df_all_annotations)}, Δt = {delta_time}s")
        video_logger.info(f"{currentFrameID}|||{delta_time:.2f}|||{len(df_current_frame)}|||{len(df_all_annotations)}")
        

    #Adding center locations
    df_all_annotations['x center'] = df_all_annotations['x top left'] + (df_all_annotations['width']/2)
    df_all_annotations['y center'] = df_all_annotations['y top left'] + (df_all_annotations['height']/2)
    
    #print(len(df_all_annotations))
    #print(df_all_annotations['width consistency'].unique())
    
    print(f'Total Objects found before video filters = {len(df_all_annotations)}')
    print(f'Saving pre-filter Video')

    # Exporting information
    root, extension = os.path.splitext(save_config['video_annotation_save_name'])
    tabular_save_name = f"{root}_pre_filtering{extension}"
    tabular_save_path = os.path.join(save_config['save_folder'], tabular_save_name) #save_config['video_annotation_save_name'])
    
    # Save file in the specified format
    _, extension = os.path.splitext(tabular_save_path)
    extension = extension[1:]
    if extension.lower() == 'parquet':
        df_all_annotations.to_parquet(tabular_save_path, index=False, engine='pyarrow')
    elif extension.lower() == 'xlsx':
        df_all_annotations.to_excel(tabular_save_path, index=False)
    else:
        raise ValueError(f"Unsupported extension '{extension}'. Use 'parquet' or 'excel'.")
    
    #Save video
    root, extension = os.path.splitext(save_config['video_save_name'])
    video_save_name = f"{root}_pre_filtering{extension}"
    video_save_path = os.path.join(save_config['save_folder'], video_save_name) #save_config['video_save_name'])
    make_annotated_video(df_all_annotations, vidcap, fps=save_config['video_fps'], 
                                           save_path=video_save_path, show_frame_number= save_config['show_frame_number'],
                                           show_object_id=save_config['show_objectID'])
    
    return df_all_annotations


def video_filters(vidcap, df, postprocess_config, search_config, save_config):
                            
    dfReturn = df.copy()
    
    # Apply percent change threshold
    """
    dfReturn['width change'] = dfReturn.groupby('objectID')['width'].pct_change()
    dfReturn['height change'] = dfReturn.groupby('objectID')['height'].pct_change()
    
    
    dfReturn = dfReturn[(dfReturn['width change'].abs() < 0.3) | (dfReturn['width change'].isna())]
    dfReturn = dfReturn[(dfReturn['height change'].abs() < 0.3) | (dfReturn['height change'].isna())]
    """
    
    #Removing short trajectories
    df_tracks_len = dfReturn.groupby('objectID')[['frame']].count()
    
    maxFrames = int(vidcap.get(cv2.CAP_PROP_FRAME_COUNT))
    min_traj_len = postprocess_config['min_traj_len']
    
    # Capping the min_traj_len
    if min_traj_len > maxFrames:
        print('min_traj_len > maxFrames... Making min_traj_len == maxFrames//4')
        min_traj_len = maxFrames//4
    
    df_tracks_len = df_tracks_len[df_tracks_len['frame']] >= min_traj_len #postprocess_config['min_traj_len']]
    objects_to_keep = df_tracks_len.index
    dfReturn = dfReturn[dfReturn['objectID'].isin(objects_to_keep)]
    
    #Removing boarder trajectories
    boarder_pixels = postprocess_config['boarder_pixels']
    firstFrame = load_frame(vidcap, frame_index=0)
    boarder_objectID = dfReturn[(dfReturn['x top left'] < boarder_pixels) |
                    (dfReturn['y top left'] < boarder_pixels) |
                    ((dfReturn['x top left'] + dfReturn['width']) > (firstFrame.shape[1] - boarder_pixels)) |
                    ((dfReturn['y top left'] + dfReturn['height']) > (firstFrame.shape[0] - boarder_pixels))
                    ]['objectID'].unique()
    dfReturn = dfReturn[~dfReturn['objectID'].isin(boarder_objectID)]
    
    #Update object id to have no missing objects
    #dfReturn = update_object_id(dfReturn)
    
    print(f'Total Objects found = {len(dfReturn)}')
    print(f'Saving post-filter Video')

    #Save tabular data
    root, extension = os.path.splitext(save_config['video_annotation_save_name'])
    tabular_save_name = f"{root}_post_filtering{extension}"
    tabular_save_path = os.path.join(save_config['save_folder'], tabular_save_name) #save_config['video_annotation_save_name'])
    
    # Save file in the specified format
    extension = extension[1:]
    if extension.lower() == 'parquet':
        dfReturn.to_parquet(tabular_save_path, index=False, engine='pyarrow')
    elif extension.lower() == 'xlsx':
        dfReturn.to_excel(tabular_save_path, index=False)
    else:
        raise ValueError(f"Unsupported extension '{extension}'. Use 'parquet' or 'excel'.")

    #Save video
    root, v_extension = os.path.splitext(save_config['video_save_name'])
    video_save_name = f"{root}_post_filtering{v_extension}"
    video_save_path = os.path.join(save_config['save_folder'], video_save_name) #save_config['video_save_name'])
    make_annotated_video(dfReturn, vidcap, fps=save_config['video_fps'], 
                                           save_path=video_save_path, show_frame_number=save_config['show_frame_number'],
                                           show_object_id=save_config['show_objectID'])
   
    #Saving additional information
    print("Computing additional information for analysis")
    dfReturn = computing_analysis_information(dfReturn, vidcap, search_config)
    
    # Save file in the specified format
    if extension.lower() == 'parquet':
        dfReturn.to_parquet(tabular_save_path, index=False, engine='pyarrow')
    elif extension.lower() == 'xlsx':
        dfReturn.to_excel(tabular_save_path, index=False)
    else:
        raise ValueError(f"Unsupported extension '{extension}'. Use 'parquet' or 'excel'.")

    print("Additional information saved")
    
    return dfReturn
    
    
def compute_stability(row, currentFrame, edge_detection_config):
    """
    Evaluate the stability of an object's bounding box by reapplying edge-based boundary detection.

    This function re-detects the boundaries of an object in the current frame using the original bounding 
    box and edge detection parameters. It computes the absolute differences between the original and 
    newly computed bounding box coordinates (top-left position and dimensions). These differences help 
    assess how stable or sensitive the object boundary is to re-detection, which can be used for 
    filtering or refining annotations.

    Parameters
    ----------
    row : pd.Series
        A row from the annotation DataFrame containing the following columns:
        - 'y top left': float
        - 'x top left': float
        - 'height': float
        - 'width': float
        - 'edgeThreshold': float
        - 'objectID': int or str

    currentFrame : np.ndarray
        The current video frame (grayscale or RGB) as a NumPy array.

    edge_detection_config : dict
        Dictionary containing edge detection parameters, specifically:
        - 'pixel_pad': int
        - optionally, 'edge_relThresh' if not using the per-row threshold

    Returns
    -------
    dict
        Dictionary with:
        - "objectID": original object ID
        - "x_diff": absolute difference in x top-left
        - "y_diff": absolute difference in y top-left
        - "w_diff": absolute difference in width
        - "h_diff": absolute difference in height
    """
    
    y_tl, x_tl, h, w, edgeThreshold, objectID = row[['y top left', 'x top left', 'height', 'width', 
                                                     'edgeThreshold', 'objectID']]
    
    y_new, x_new, h_new, w_new, madeBy, edgeThreshold, nms = detect_and_refine_boundary(
        currentFrame, x_tl, y_tl, w, h,
        pixel_pad=edge_detection_config['pixel_pad'],
        edge_relThresh= edgeThreshold, #edge_detection_config['edge_relThresh'],
    )

    return {
        "objectID": objectID,
        "x_diff": abs(x_tl - x_new),
        "y_diff": abs(y_tl - y_new),
        "w_diff": abs(w - w_new),
        "h_diff": abs(h - h_new)
    }


def remove_unstable(currentFrame, df, edge_detection_config, pixel_threshold=None, parallel_processing=False):
    if pixel_threshold is None:
        pixel_threshold = math.ceil(edge_detection_config['pixel_pad']/2)
    
    if parallel_processing and len(df) > 1:
        jobs = [
            delayed(compute_stability)(row, currentFrame, edge_detection_config)
            for _, row in df.iterrows()
        ]
        refined_rows = run_balanced_parallel(jobs)
    else:
        refined_rows = [
            compute_stability(row, currentFrame, edge_detection_config)
            for _, row in df.iterrows()
        ]
    
    # Step 2: Filter stable boxes
    refined_df = pd.DataFrame(refined_rows)
    stable_ids = refined_df[
        (refined_df["x_diff"] < pixel_threshold) &
        (refined_df["y_diff"] < pixel_threshold) &
        (refined_df["w_diff"] < pixel_threshold) &
        (refined_df["h_diff"] < pixel_threshold)
    ]["objectID"].tolist()

    # Step 3: Return only stable rows    
    return df[df["objectID"].isin(stable_ids)].reset_index(drop=True)
    
    
def computing_analysis_information(df, vidcap, search_config):
    df = initial_annotation_info(df)
    df = compute_crop_information(df, vidcap, search_config)
    df = compute_displacement(df)

    df = df.sort_values(by=['frame', 'objectID'])
    
    return df

def compute_displacement(df):
    # Sort to ensure correct previous-frame access
    df = df.sort_values(by=['objectID', 'frame'])

    # Compute x and y displacements using groupby + diff
    df['x_displacement'] = df.groupby('objectID')['x top left'].diff()
    df['y_displacement'] = df.groupby('objectID')['y top left'].diff()

    # Optional: Compute Euclidean distance moved (optional)
    df['euclidean_distance'] = (df['x_displacement']**2 + df['y_displacement']**2)**0.5

    return df


def initial_annotation_info(df):
    # Add information about initial annotation
    first_frame_human_flags = (df[df['frame'] == 0].set_index('objectID')['humanMade']
                               .rename('Initial Annotation Madeby')
                              )
    df['Initial Annotation Madeby'] = (df['objectID'].map(first_frame_human_flags)
                                              .replace({True, False}, {'Manual', 'Algorithmic'}))

    return df
    
def compute_crop_features(row, currentFrame):
    y1 = math.floor(row["y top left"])
    y2 = math.ceil(row["y top left"] + row["height"])
    x1 = math.floor(row["x top left"])
    x2 = math.ceil(row["x top left"] + row["width"])
    
    im_crop = crop(currentFrame, y1, y2, x1, x2)
    h_com, w_com = center_of_mass(normalizeImage(im_crop))

    centroid_y = (x2 - x1) / 2
    centroid_x = (y2 - y1) / 2

    return {
        'Centroid Y': centroid_y,
        'Centroid X': centroid_x,
        'Center of mass Y': h_com,
        'Center of mass X': w_com,
        'entropy': shannon_entropy(im_crop),
        'Peak Intensity': np.max(im_crop),
        'Mean Intensity': np.mean(im_crop),
        'Local Contrast': np.std(im_crop)
    }	
	
def compute_crop_information(df, vidcap, search_config):
    frames = sorted(df['frame'].unique())
    
    for frame_value in frames:
        mask = df['frame'] == frame_value
        df_temp = df[mask].copy()
        currentFrame = load_frame(vidcap, frame_index=frame_value)
        
        if search_config['parallel_processing']:
            jobs = [
                delayed(compute_crop_features)(row, currentFrame)
                for _, row in df_temp.iterrows()
            ]
            results = run_balanced_parallel(jobs)
        else:
            results = [
                compute_crop_features(row, currentFrame)
                for _, row in df_temp.iterrows()
            ]
			
        results_df = pd.DataFrame(results, index=df_temp.index)
        df.loc[mask, results_df.columns] = results_df
    
    # Post-processing
    df['Center Uncertainty Y'] = abs(df['Center of mass Y'] - df['Centroid Y'])
    df['Center Uncertainty X'] = abs(df['Center of mass X'] - df['Centroid X'])
    df['Total Uncertainty'] = df['Center Uncertainty Y'] + df['Center Uncertainty X']
    df['Percent Uncertainty'] = 100 * (
        (df['Center Uncertainty Y'] / df['height']) +
        (df['Center Uncertainty X'] / df['width'])
    )
    
    return df


def create_mask_from_annotations(df: pd.DataFrame, image_shape: tuple, pad: int = 0) -> np.ndarray:
    """
    Create a binary mask from object annotations.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing 'x top left', 'y top left', 'width', and 'height'.
    image_shape : tuple
        Shape of the original image (H, W).
    pad : int
        Padding to expand the bounding boxes.
    
    Returns:
    --------
    mask : np.ndarray
        Binary mask with ones where objects are located.
    """
    mask = np.zeros(image_shape[:2], dtype=np.uint8)
    
    for _, row in df.iterrows():
        x1 = max(int(row['x top left']) - pad, 0)
        y1 = max(int(row['y top left']) - pad, 0)
        x2 = min(int(row['x top left'] + row['width']) + pad, image_shape[1])
        y2 = min(int(row['y top left'] + row['height']) + pad, image_shape[0])
        
        mask[y1:y2, x1:x2] = 1

    return mask
    
    
def apply_mask_to_image(image: np.ndarray, mask: np.ndarray, mode='black') -> np.ndarray:
    """
    Apply a binary mask to an image.

    Parameters:
    -----------
    image : np.ndarray
        Original image.
    mask : np.ndarray
        Binary mask.
    mode : str
        If 'black', zero out the masked areas. If 'blur', apply a blur.

    Returns:
    --------
    modified_image : np.ndarray
        Masked image.
    """
    modified_image = image.copy()
    
    if mode.lower() == 'black':
        modified_image[mask == 1] = 0
    elif mode.lower() == 'none':
        modified_image[mask == 1] = 1
    elif mode.lower() == 'white':
        modified_image[mask == 1] = 255
    elif mode.lower() == 'blur':
        blurred = cv2.GaussianBlur(image, (15, 15), 0)
        modified_image[mask == 1] = blurred[mask == 1]
    else:
        raise ValueError(f"Unsupported mode '{extension}'. Use 'black', 'white', 'blur', or 'none'.")

    return modified_image

