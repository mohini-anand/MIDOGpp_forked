def toy_locations(particle_shape = "Gaussian", box = "Overfit", objectNum=0):
    madeBy = 'Human'
    if particle_shape == "Gaussian" and box == "Overfit" and objectNum==0:
        #over sized box, Gaussian: A bit imperfect dimension
        x_tl, y_tl, w, h = 187, 86, 16, 17
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Gaussian" and box == "Overfit" and objectNum==1:
        #over sized box, Gaussian: A bit imperfect dimension
        x_tl, y_tl, w, h = 136, 486, 12, 12 #378, 378, 13, 13
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Gaussian" and box == "Overfit" and objectNum==2:
        #over sized box, Gaussian: A bit imperfect dimension
        x_tl, y_tl, w, h = 250, 209, 8, 8
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Gaussian" and box == "Overfit" and objectNum==3:
        #Not true Gaussian
        x_tl, y_tl, w, h = 76, 90, 25, 25
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Gaussian" and box == "Overfit" and objectNum==4:
        #Not true Gaussian
        x_tl, y_tl, w, h = 140, 220, 15, 15
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
        
    elif particle_shape == "Defocused" and box == "Overfit" and objectNum==0:
        #over sized box, Defocused
        x_tl, y_tl, w, h = 233, 336, 262-233, 363-336 #severe
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Gaussian" and box == "Underfit" and objectNum==0:
        #under sized box, Gaussian
        #x_tl, y_tl, w, h = 187, 86, 16, 9
        x_tl, y_tl, w, h = 187, 86, 16, 4
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Defocused" and box == "Underfit" and objectNum==0:
        #under sized box, Defocused
        #x_tl, y_tl, w, h = 233, 336, 262-233, 363-346
        x_tl, y_tl, w, h = 233, 336, 262-233, 363-356 #severe
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Gaussian" and box == "Corner" and objectNum==0:
        #Particle on bottom left corner
        x_tl, y_tl, w, h = 192, 80, 12, 15
        #x_tl, y_tl, w, h = 191, 80, 12, 15
        print(f"Original dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Defocused" and box == "Corner" and objectNum==0:
        #under sized box, Defocused
        x_tl, y_tl, w, h = 242, 330, 25, 15
        #x_tl, y_tl, w, h = 240, 330, 25, 15
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Defocused" and box == "Overlapping" and objectNum==0:
        #under sized box, Defocused
        x_tl, y_tl, w, h = 90, 320, 70, 70
        #x_tl, y_tl, w, h = 240, 330, 25, 15
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    elif particle_shape == "Defocused" and box == "Overlapping" and objectNum==1:
        #under sized box, Defocused
        x_tl, y_tl, w, h = 62, 64, 15, 25
        #x_tl, y_tl, w, h = 240, 330, 25, 15
        print(f"Initial dimensions: {x_tl = }, {y_tl = }, {w = }, {h = }")
    else:
        print('Wrong selection of particle_shape or box')
        return 0, 0, 0, 0
    
    return y_tl, x_tl, h, w, madeBy #x_tl, y_tl, w, h, madeBy
    