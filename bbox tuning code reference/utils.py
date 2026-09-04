import os
import csv

import pandas as pd
import numpy as np
import cv2

import time
import math

import matplotlib as mpl
from matplotlib import pyplot as plt
import plotly.express as px

import gc
import logging

from skimage.filters import laplace
from skimage.measure import shannon_entropy


def plot_tps(image_x, image_y, image_threshold): 
    
    plt.figure(figsize=(12,8))
    plt.plot(list(range(len(image_x))), image_x, color ="red") 
    plt.plot(list(range(len(image_y))), image_y, color ="green") 
    
    plt.legend(["x direction", "y direction"])
    
    x1_loc =  np.argwhere(image_x>image_threshold)[0][0]
    x2_loc =  np.argwhere(image_x>image_threshold)[-1][0]
    y1_loc =  np.argwhere(image_y>image_threshold)[0][0]
    y2_loc =  np.argwhere(image_y>image_threshold)[-1][0]
    
    plt.plot(x1_loc, image_x[x1_loc], markersize=18, color = "red", 
             marker=mpl.markers.CARETRIGHTBASE) # x_tl+w
    plt.plot(y1_loc, image_y[y1_loc], markersize=18, color = "green", 
             marker=mpl.markers.CARETRIGHTBASE) # x_tl+w
    plt.plot(x2_loc, image_x[x2_loc], markersize=18, color = "red", 
             marker=mpl.markers.CARETLEFTBASE) # x_tl+w
    plt.plot(y2_loc, image_y[y2_loc], markersize=18, color = "green", 
             marker=mpl.markers.CARETLEFTBASE) # x_tl+w
    
    plt.show()
    
    return None


def update_object_id(df: pd.DataFrame) -> pd.DataFrame:
    """
    Renumber object IDs in the DataFrame to be consecutive integers starting from 0,
    preserving the original order of appearance.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with a column named 'objectID'.

    Returns:
    --------
    pd.DataFrame
        DataFrame with updated 'objectID' values.
    """
    # Get unique sorted values from the column
    unique_sorted_values = sorted(df['objectID'].unique())

    # Create a mapping from old IDs to new consecutive ones
    value_mapping = {val: idx for idx, val in enumerate(unique_sorted_values)}

    # Apply the mapping to the 'objectID' column
    df['objectID'] = df['objectID'].map(value_mapping)

    return df


def make_annotated_video(
    df: pd.DataFrame,
    vidcap,
    frame_skip: int = 1,
    fps: int = 15,
    save_path: str = 'output.mp4',
    show_frame_number: bool = False,
    show_object_id: bool = False
) -> np.ndarray:
    """
    Generate an annotated video from a video capture and bounding box dataframe.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with bounding box annotations. Must include 'frame', 'x top left', 'y top left', 'width', 'height'.
    vidcap : cv2.VideoCapture
        OpenCV video capture object.
    frame_skip : int
        Process every nth frame.
    fps : int
        Frames per second for output video.
    save_path : str
        File path to save the output video.
    show_frame_number : bool
        Overlay the frame index on the top-right corner.
    show_object_id : bool
        Overlay object ID above each bounding box.

    Returns:
    --------
    np.ndarray
        The last frame written to the video (for optional preview).
    """
    frames_to_visualize = sorted([f for f in df['frame'].unique() if f % frame_skip == 0])
    if not frames_to_visualize:
        raise ValueError("No frames selected for visualization.")

    # Determine video resolution from the first valid frame
    first_frame = get_frame(vidcap, frame_index=frames_to_visualize[0])
    height, width = first_frame.shape[:2]
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(save_path, fourcc, fps, (width, height))
    
    font = cv2.FONT_HERSHEY_SIMPLEX
    
    # Rename columns for dot-access
    df = df.rename(columns={"x top left": "x", "y top left": "y",
        "width": "w", "height": "h"})
    
    #Vectorized Grouping according to the frames (avoid re-filtering in loop)
    frame_groups = dict(tuple(df.groupby("frame")))
    
    for frame_index in frames_to_visualize:
        current_frame = load_frame(vidcap, frame_index=frame_index)
        #df_temp = df[df['frame'] == frame_index]
        df_temp = frame_groups.get(frame_index, pd.DataFrame())
        img_color = cv2.cvtColor(current_frame, cv2.COLOR_GRAY2RGB)

        # Draw bounding boxes and IDs
        for row in df_temp.itertuples():
            #x, y, w, h = map(int, [row["x top left"], row["y top left"], row["width"], row["height"]])
            x, y, w, h = int(row.x), int(row.y), int(row.w), int(row.h)
            top_left = (x, y)
            bottom_right = (x + w, y + h)
            cv2.rectangle(img_color, top_left, bottom_right, (0, 0, 255), 1)

            if show_object_id:
                text_x = max(0, int(x))
                text_y = max(0, int(y - 5))
                cv2.putText(
                    img_color,
                    str(row.objectID),
                    (text_x, text_y),
                    font,
                    0.4,
                    (0, 0, 255),
                    1,
                    cv2.LINE_AA,
                )

        # Add frame number in top-right
        if show_frame_number:
            text = f"{int(frame_index)}"
            font = cv2.FONT_HERSHEY_SIMPLEX
            font_scale = 1
            font_thickness = 2
            text_size, _ = cv2.getTextSize(text, font, font_scale, font_thickness)
            text_x = width - text_size[0] - 5
            text_y = 30
            cv2.putText(img_color, text, (text_x, text_y), font, font_scale, (0, 0, 255), font_thickness, cv2.LINE_AA)

        out.write(img_color)

    out.release()
    cv2.destroyAllWindows()
    gc.collect()

    print(f"Video saved to {save_path}")
    return None


def center_of_mass(arr: np.ndarray) -> tuple[float, float]:
    """
    Compute the center of mass (CoM) of a 2D array, where intensity values represent "mass".

    Parameters:
    -----------
    arr : np.ndarray
        2D array of intensity values (e.g., grayscale image or binary mask).

    Returns:
    --------
    (height_com, width_com) : tuple[float, float]
        Vertical and horizontal coordinates of the center of mass.
    """
    total_mass = np.sum(arr)
    
    # Calculate vertical (row-wise) and horizontal (col-wise) center of mass
    height_indices = np.arange(arr.shape[0])[:, np.newaxis]  # column vector
    width_indices = np.arange(arr.shape[1])[np.newaxis, :]   # row vector

    height_com = np.sum(arr * height_indices) / total_mass
    width_com = np.sum(arr * width_indices) / total_mass

    return height_com, width_com


def get_locations(
    y_tl: float,
    x_tl: float,
    h: float,
    w: float,
    box_shake_values: list,
    pixels_to_move: int
) -> pd.DataFrame:
    """
    Generate jittered versions of a bounding box by shifting its position.

    Parameters:
    -----------
    y_tl : float
        Top-left y-coordinate of the original bounding box.
    x_tl : float
        Top-left x-coordinate of the original bounding box.
    h : float
        Height of the bounding box.
    w : float
        Width of the bounding box.
    box_shake_values : list of int
        Multipliers to define how far to jitter the box.
    pixels_to_move : int
        Base number of pixels used for shifting in each direction.

    Returns:
    --------
    pd.DataFrame
        DataFrame of jittered bounding box coordinates.
        Columns: ['y_tl', 'x_tl', 'h', 'w']
    """
    if box_shake_values is None:
        box_shake_values = [0]

    y_tl, x_tl, h, w = float(y_tl), float(x_tl), float(h), float(w)
    locations = [{'y_tl': y_tl, 'x_tl': x_tl, 'h': h, 'w': w}]  # Original box

    # Add jittered variants
    for multiplier in box_shake_values:
        offset = pixels_to_move * multiplier
        shifts = [
            ( y_tl - offset, x_tl             ),
            ( y_tl + offset, x_tl             ),
            ( y_tl,           x_tl - offset   ),
            ( y_tl,           x_tl + offset   ),
            ( y_tl - offset, x_tl - offset    ),
            ( y_tl - offset, x_tl + offset    ),
            ( y_tl + offset, x_tl - offset    ),
            ( y_tl + offset, x_tl + offset    ),
        ]
        for y, x in shifts:
            locations.append({'y_tl': y, 'x_tl': x, 'h': h, 'w': w})

    return pd.DataFrame(locations).drop_duplicates().reset_index(drop=True)


#Might be a redundant function (or have redundant features)
def visualize_frames(
    df: pd.DataFrame,
    vidcap,
    frame_skip: int = 10
) -> None:
    """
    Visualize annotated frames from a video at regular intervals.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing bounding box information.
        Required columns: 'frame', 'x top left', 'y top left', 'width', 'height'.
    vidcap : cv2.VideoCapture
        OpenCV video capture object.
    frame_skip : int
        Number of frames to skip between visualizations.

    Returns:
    --------
    None
    """
    frames_to_visualize = [
        frame for frame in df['frame'].unique()
        if frame % frame_skip == 0
    ]

    for frame_index in sorted(frames_to_visualize):
        # Load the frame and get its annotations
        frame = load_frame(vidcap, frame_index=frame_index)
        frame_rgb = cv2.cvtColor(frame, cv2.COLOR_GRAY2RGB)
        df_frame = df[df['frame'] == frame_index]

        # Plot with matplotlib
        fig, ax = plt.subplots(figsize=(15, 15))
        ax.imshow(frame_rgb)

        for _, row in df_frame.iterrows():
            rect = mpl.patches.Rectangle(
                (math.floor(row["x top left"]), math.floor(row["y top left"])),
                math.ceil(row["width"]),
                math.ceil(row["height"]),
                linewidth=1,
                edgecolor='r',
                facecolor='none'
            )
            ax.add_patch(rect)

        ax.set_title(f"Frame {frame_index}")
        plt.show()

    
def annotate_bounding_boxes_over_image(
    image: np.ndarray,
    df: pd.DataFrame,
    add_object_id: bool = False
) -> np.ndarray:
    """
    Draw bounding boxes on an image using coordinates from a DataFrame.

    Parameters:
    -----------
    image : np.ndarray
        Input grayscale or RGB image.
    df : pd.DataFrame
        DataFrame with columns: 'x top left', 'y top left', 'width', 'height', and optionally 'objectID'.
    add_object_id : bool
        Whether to annotate each box with its objectID.

    Returns:
    --------
    np.ndarray
        Image with bounding boxes drawn on it.
    """
    # Convert grayscale to RGB for colored annotations
    if len(image.shape) == 2:
        image = cv2.cvtColor(image, cv2.COLOR_GRAY2RGB)

    for _, row in df.iterrows():
        x1 = int(math.floor(row["x top left"]))
        y1 = int(math.floor(row["y top left"]))
        x2 = int(math.ceil(row["x top left"] + row["width"]))
        y2 = int(math.ceil(row["y top left"] + row["height"]))

        # Draw rectangle
        cv2.rectangle(image, (x1, y1), (x2, y2), color=(0, 0, 255), thickness=2)

        # Draw object ID if requested
        if add_object_id and "objectID" in row:
            cv2.putText(
                image,
                str(row["objectID"]),
                (x1, max(0, y1 - 5)),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.8,
                (0, 0, 255),
                2,
                cv2.LINE_AA,
            )

    return image


def mask(
    image: np.ndarray,
    row: pd.Series,
    offset: int = 0
) -> np.ndarray:
    """
    Apply a rectangular binary mask around an object defined in a DataFrame row.

    Parameters:
    -----------
    image : np.ndarray
        Input image (grayscale or RGB).
    row : pd.Series
        Row containing 'x top left', 'y top left', 'width', 'height' fields.
    offset : int
        Number of pixels to expand the mask region on all sides.

    Returns:
    --------
    np.ndarray
        Masked image with values preserved only in the ROI.
    """
    height, width = image.shape[:2]

    # Compute crop bounds with clipping
    top = max(0, int(row["y top left"] - offset))
    bottom = min(height, int(row["y top left"] + row["height"] + offset))
    left = max(0, int(row["x top left"] - offset))
    right = min(width, int(row["x top left"] + row["width"] + offset))

    # Initialize mask and apply region
    mask = np.zeros(image.shape[:2], dtype=np.uint8)
    mask[top:bottom, left:right] = 1

    # Apply mask (preserves original image dtype)
    return cv2.bitwise_and(image, image, mask=mask)

    
def load_frame(
    vidcap: cv2.VideoCapture,
    frame_index: int = 0,
    display_frame: bool = False
) -> np.ndarray:
    """
    Load a preprocessed grayscale (green channel) frame from the video.

    Parameters:
    -----------
    vidcap : cv2.VideoCapture
        OpenCV video capture object.
    frame_index : int
        Index of the frame to load.
    display_frame : bool
        Whether to display the raw frame before processing.

    Returns:
    --------
    np.ndarray
        Preprocessed grayscale image (green channel only).
    """
    # Retrieve the frame using get_frame
    frame = get_frame(vidcap, frame_index=frame_index, display_frame=display_frame)
    if frame is None:
        raise ValueError(f"Could not load frame at index {frame_index}.")

    # Extract green channel only
    green_channel = frame[:, :, 1]

    # Apply preprocessing
    processed = pre_processing(green_channel)

    return processed


def pre_processing(image: np.ndarray, kernel_size = 61) -> np.ndarray:
    """
    Apply background correction, local contrast enhancement, and sharpening to a grayscale image.

    This preprocessing pipeline is useful for preparing microscopy or low-contrast grayscale images 
    by removing background illumination bias, enhancing local contrast, and sharpening edges. 

    It performs the following steps:
    1. Background estimation and removal using morphological opening.
    2. Local contrast enhancement by subtracting background from the original image.
    3. Gentle sharpening using unsharp masking with Gaussian blur.

    Parameters
    ----------
    image : np.ndarray
        Input image in 8-bit, single-channel grayscale format (dtype=np.uint8).

    kernel_size : int, optional
        Diameter of the elliptical kernel used for morphological opening. 
        A larger value suppresses larger background features (default is 61). 
        Must be a positive odd integer.

    Returns
    -------
    np.ndarray
        Output image after preprocessing, with reduced background, enhanced contrast, and sharpened features.
    """
    # Step 1: Estimate the background using morphological opening
    background = cv2.morphologyEx(image, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size)))

    # Step 2: Subtract background
    corrected = cv2.subtract(image, background)

    # Optional: Gentle denoise and sharpen (your previous logic)
    blurred = cv2.GaussianBlur(corrected, (11, 11), sigmaX=2.0)
    sharpened = cv2.addWeighted(corrected, 1.5, blurred, -0.5, 0)
    
    return sharpened
    

def get_frame(vidcap, frame_index=0, display_frame=False):
    """
    Extract a specific frame from a video capture object.

    Parameters:
    -----------
    vidcap : cv2.VideoCapture
        OpenCV video capture object.
    frame_index : int
        Index of the frame to extract.
    display_frame : bool
        Whether to display the extracted frame using display_image.

    Returns:
    --------
    frame : np.ndarray or None
        Extracted frame as an image (BGR format). Returns None if frame_index is invalid.
    """
    total_frames = int(vidcap.get(cv2.CAP_PROP_FRAME_COUNT))
    frame_index = int(frame_index)

    if 0 <= frame_index < total_frames:
        # Set the current frame position
        vidcap.set(cv2.CAP_PROP_POS_FRAMES, frame_index)
    else:
        print(f"[ERROR] Frame index must be between 0 and {total_frames - 1}, but got {frame_index}.")
        return None

    ret, frame = vidcap.read()
    
    if not ret:
        print("[ERROR] Failed to read the frame.")
        return None

    if display_frame:
        print(f"[INFO] Displaying frame {frame_index}")
        print("Frame shape:", frame.shape)
        display_image(frame)

    # Reset video to the beginning
    vidcap.set(cv2.CAP_PROP_POS_FRAMES, 0)

    return frame


def uncertainty_below_limit(df_use, max_center_uncertainty):
    """
    Returns a summary string indicating how many objects have 
    center uncertainty below a given threshold in both X and Y directions.

    Parameters:
    -----------
    df_use : pd.DataFrame
        DataFrame containing columns 'Center Uncertainty X' and 'Center Uncertainty Y'.
    max_center_uncertainty : float
        Threshold for maximum acceptable uncertainty.

    Returns:
    --------
    str
        Summary text indicating number and percentage of valid objects.
    """
    # Filter the DataFrame for rows that meet both X and Y uncertainty thresholds
    mask = (
        (df_use['Center Uncertainty Y'] <= max_center_uncertainty) & 
        (df_use['Center Uncertainty X'] <= max_center_uncertainty)
    )
    df_temp = df_use[mask]

    text = '{}/{} ({}%) of objects have max_center_uncertainty <= {} in both directions'.format(
        len(df_temp), len(df_use), round(100 * len(df_temp) / len(df_use), 2), max_center_uncertainty
    )

    return text


def crop(img, vCropStart, vCropEnd, hCropStart, hCropEnd):
    """
    Safely crop a region from the image, ensuring indices remain within bounds.

    Parameters:
    -----------
    img : np.ndarray
        Input image (grayscale or RGB).
    vCropStart, vCropEnd : int
        Vertical crop start and end indices (rows).
    hCropStart, hCropEnd : int
        Horizontal crop start and end indices (columns).

    Returns:
    --------
    np.ndarray
        Cropped region of the image, guaranteed to be valid.
    """
    # Get image dimensions
    height, width = img.shape[:2]

    # Ensure vertical crop start is at least 0 and no more than the second-last row
    vCropStart = max(0, min(vCropStart, height - 2))

    # Ensure vertical crop end is at least one row after start and within image bounds
    vCropEnd = max(vCropStart + 1, min(vCropEnd, height - 1))

    # Ensure horizontal crop start is at least 0 and no more than the second-last column
    hCropStart = max(0, min(hCropStart, width - 2))

    # Ensure horizontal crop end is at least one column after start and within image bounds
    hCropEnd = max(hCropStart + 1, min(hCropEnd, width - 1))

    # Return the cropped image slice
    return img[vCropStart:vCropEnd, hCropStart:hCropEnd]


def row_crop(
    image: np.ndarray,
    row: pd.Series,
    offset: int = 0
) -> np.ndarray:
    """
    Crop a region from the image using bounding box information from a DataFrame row.

    Parameters:
    -----------
    image : np.ndarray
        Input image (grayscale or color).
    row : pd.Series
        A row containing bounding box information with keys:
        'x top left', 'y top left', 'width', 'height'.
    offset : int
        Optional padding (in pixels) to expand the crop on all sides.

    Returns:
    --------
    np.ndarray
        Cropped image region.
    """
    try:
        x_tl = row["x top left"]
        y_tl = row["y top left"]
        w = row["width"]
        h = row["height"]
    except KeyError as e:
        raise KeyError(f"Missing bounding box key: {e}")

    top = math.floor(y_tl - offset)
    bottom = math.ceil(y_tl + h + offset)
    left = math.floor(x_tl - offset)
    right = math.ceil(x_tl + w + offset)

    return crop(image, top, bottom, left, right)


def padded_crop(
    image_padded: np.ndarray,
    vCropStart: int,
    vCropEnd: int,
    hCropStart: int,
    hCropEnd: int,
) -> np.ndarray:
    """
    Crop a padded image region.

    Parameters:
    -----------
    image_padded : np.ndarray
        Input image with padding already applied.
    vCropStart : int
        Starting vertical pixel index.
    vCropEnd : int
        Ending vertical pixel index.
    hCropStart : int
        Starting horizontal pixel index.
    hCropEnd : int
        Ending horizontal pixel index.

    Returns:
    --------
    np.ndarray
        Cropped image region.
    """
    return crop(image_padded, vCropStart, vCropEnd, hCropStart, hCropEnd)


def display_image(
    image: np.ndarray,
    axis: bool = True,
    title: str = "",
    style: str = "mpl",
    figsize: tuple = (15, 15),
    save_path: str = None
) -> None:
    """
    Display a grayscale image using matplotlib or plotly, and optionally save it.

    Parameters:
    -----------
    image : np.ndarray
        Input grayscale or multi-channel image.
    axis : bool
        Whether to show axis ticks (only applies to matplotlib).
    title : str
        Title to be displayed on the image.
    style : str
        Display method: 'mpl' (matplotlib) or 'plotly'.
    figsize : tuple
        Size of the figure (for matplotlib).
    save_path : str
        File path to save the image. Saves in matplotlib style regardless of display mode.
    """
    print("Image Shape:", image.shape)

    # Save image using matplotlib backend, regardless of style
    if save_path:
        plt.figure(figsize=figsize)
        if not axis:
            plt.axis('off')
        plt.title(title)
        if len(image.shape) == 2:
            plt.imshow(image, cmap='gray')
        else:
            plt.imshow(image)
        plt.tight_layout()
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

    # Display logic
    if style == "mpl":
        plt.figure(figsize=figsize)
        if not axis:
            plt.axis('off')
        plt.title(title)
        if len(image.shape) == 2:
            plt.imshow(image, cmap='gray')
        else:
            plt.imshow(image)
        plt.show()
        
    elif style == "plotly":
        if len(image.shape) == 2:
            fig = px.imshow(image, title=title, color_continuous_scale="gray")
        else:
            fig = px.imshow(image, title=title)
        fig.update_layout(width=figsize[0]*60, height=figsize[1]*60)
        fig.show()

    else:
        raise ValueError(f"Unknown style: {style}. Use 'mpl' or 'plotly'.")


def display_crop(
    image: np.ndarray,
    y_tl: float,
    x_tl: float,
    h: float,
    w: float,
    style: str = "mpl"
) -> None:
    """
    Crop a region from the image and display it.

    Parameters:
    -----------
    image : np.ndarray
        Input image (grayscale or RGB).
    y_tl, x_tl : float
        Top-left coordinates of the bounding box.
    h, w : float
        Height and width of the bounding box.
    style : str
        Display method: 'mpl' for matplotlib, 'plotly' for interactive.
    """
    im_crop = crop(
        image,
        math.floor(y_tl), math.ceil(y_tl + h),
        math.floor(x_tl), math.ceil(x_tl + w)
    )

    if style == "mpl":
        display_image(im_crop)
    elif style == "plotly":
        fig = px.imshow(im_crop)
        fig.show()
    else:
        raise ValueError(f"Unknown display style: '{style}'")

    return None


def display_annotations_on_frame(
    image: np.ndarray,
    df: pd.DataFrame,
    display_number: bool = False,
    save_image: bool = False,
    save_name: str = None,
    save_folder: str = "output_files"
) -> None:
    """
    Display bounding boxes on an image using matplotlib.

    Parameters:
    ----------
    image : np.ndarray
        Grayscale or RGB image array.

    df : pd.DataFrame
        DataFrame containing bounding box data.
        Required columns: 'x top left', 'y top left', 'width', 'height', 'humanMade', 'objectID', 'frame'.

    display_number : bool, default=False
        Whether to overlay object IDs on the image.

    save_image : bool, default=False
        Whether to save the annotated image.

    save_folder : str, default="output_files"
        Directory to save the annotated image if save_image is True.
    """
    if df.empty:
        print("Warning: DataFrame is empty. Nothing to annotate.")
        return

    # Ensure image is 3-channel RGB
    if len(image.shape) == 2:
        image = cv2.cvtColor(image, cv2.COLOR_GRAY2RGB)

    fig, ax = plt.subplots(figsize=(15, 15))
    ax.imshow(image)

    for _, row in df.iterrows():
        x, y = math.floor(row["x top left"]), math.floor(row["y top left"])
        w, h = math.ceil(row["width"]), math.ceil(row["height"])

        # Color and linestyle by annotation type
        is_human = row.get("humanMade", False)
        color = 'y' if is_human else 'r'
        linestyle = '-' if is_human else '--'

        rect = mpl.patches.Rectangle((x, y), w, h,
                                     linewidth=1,
                                     edgecolor=color,
                                     linestyle=linestyle,
                                     facecolor='none')
        ax.add_patch(rect)

        if display_number and "objectID" in row:
            ax.text(x, y - 5, str(row["objectID"]),
                    fontsize=10, color=color, fontweight='bold')

    if save_image:
        os.makedirs(save_folder, exist_ok=True)
        frame_id = df['frame'].iloc[0] if 'frame' in df.columns else 'unknown'
        if save_name is None:
            save_name = f"frame{frame_id}_annotations{len(df)}_{x}.jpg"
        save_path = os.path.join(save_folder, save_name)

        plt.axis('off')
        fig.savefig(save_path, bbox_inches='tight', dpi=450, pad_inches=0)
        print(f"Saved annotated image to: {save_path}")
        plt.axis('on')

    plt.show()



def display_individual_crops(df: pd.DataFrame, image: np.ndarray, style: str = 'mpl', 
    figsize: tuple = (3, 3)) -> None:
    """
    Display individual cropped images of objects/objects from a given frame.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing bounding box information with at least:
        'x top left', 'y top left', 'width', 'height', and 'objectID'.

    currentFrame : np.ndarray
        Full frame (image) from which crops will be extracted.

    style : str, default='mpl'
        Display style passed to `display_image`. Options may include 'mpl', 'plotly', etc.

    figsize : tuple, default=(3, 3)
        Size of the display figure for each crop.
    """
    if df.empty:
        print("Warning: The DataFrame is empty. No crops to display.")
        return

    for _, row in df.iterrows():
        print(f"Displaying object: {row['objectID']}")
        crop = row_crop(image, row, offset=0)
        display_image(crop, style=style, figsize=figsize)

    return None


def normalizeImage(image: np.ndarray) -> np.ndarray:
    """
    Normalize an image to 8-bit grayscale (0–255) range.

    Parameters:
    -----------
    image : np.ndarray
        Input image (grayscale or single-channel float/int image).

    Returns:
    --------
    np.ndarray
        Normalized 8-bit image.
    """
    if image is None:
        raise ValueError("Input image cannot be None.")
    
    #Avoiding unexpected truncation from np.uint8
    image = image.astype(np.float32)
    min_val = image.min()
    max_val = image.max()

    # Prevent division by zero
    if max_val == min_val:
        return np.zeros_like(image, dtype=np.uint8)

    norm = (image - min_val) / (max_val - min_val)
    
    return (norm * 255).astype(np.uint8)
    
    
class CSVMetricLogger(logging.Handler):
    def __init__(self, csv_path, header, delimiter="|||", datefmt="%Y-%m-%d %H:%M:%S"):
        """
        Parameters:
        - csv_path: path to output CSV
        - header: list of column names. First two must be timestamp and log level
        - delimiter: log field separator (default '|||')
        - datefmt: timestamp format
        """
        super().__init__()
        self.csv_path = csv_path
        self.header = header
        self.delimiter = delimiter
        self.datefmt = datefmt

        # Create CSV file with header if it doesn't exist
        if not os.path.exists(self.csv_path):
            with open(self.csv_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(self.header)

    def emit(self, record):
        try:
            # Split the message body using delimiter
            message_parts = record.getMessage().split(self.delimiter)
            #print(f'{message_parts = }')
            
             # Format timestamp using Formatter instance
            if not hasattr(self, '_formatter'):
                self._formatter = logging.Formatter(datefmt=self.datefmt)
            timestamp = self._formatter.formatTime(record)

            # Create full row: timestamp + level + message parts
            row = [timestamp, record.levelname]
            #print(f"{row = }")
            row.extend(part.strip() for part in message_parts)

            # Pad row if message is short
            while len(row) < len(self.header):
                row.append("")

            # Truncate if message is too long
            row = row[:len(self.header)]

            # Write to CSV
            with open(self.csv_path, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(row)

        except Exception:
            self.handleError(record)
            

def setup_logger(name='unnamed_log', video_name='', log_folder='log_files', 
                 log_header="Log Time|||Log Level|||Log Location|||Time (s)|||Rows found\n",
                 datefmt='%Y-%m-%d %H:%M:%S'):
    """
    Set up a logger with both standard and CSV file handlers.

    Parameters:
    - name (str): Logger name.
    - video_name (str): Name of the video file (used to name the log files).
    - log_folder (str): Folder where logs will be saved.
    - log_header (str): Header line for log files, separated by '|||'.
    - datefmt (str): Date format used in timestamps.

    Returns:
    - logging.Logger: Configured logger instance.
    """
    # Ensure log folder exists
    os.makedirs(log_folder, exist_ok=True)
    
    # Log file paths
    base_name = video_name.split('.')[0]
    log_path = os.path.join(log_folder, f'{base_name}_{name}.log')
    csv_log_path = os.path.join(log_folder, f'{base_name}_{name}_log.csv')
    csv_header = log_header.strip().split('|||')

    # Prepare logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()

    # Write headers
    try:
        with open(log_path, "w") as f:
            f.write(log_header)
    except IOError as e:
        print(f"Failed to write log header to {log_path}: {e}")

    # Formatter
    formatter = logging.Formatter(fmt='%(asctime)s|||%(levelname)s|||%(message)s', datefmt=datefmt)

    # File handler
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)

    # CSV handler (assumes CSVMetricLogger is defined elsewhere)
    csv_handler = CSVMetricLogger(csv_log_path, header=csv_header, datefmt=datefmt)

    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(csv_handler)

    return logger