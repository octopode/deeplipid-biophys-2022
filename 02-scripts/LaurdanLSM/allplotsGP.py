import numpy as np
import pandas as pd
import seaborn as sns
from skimage import data, color, io
import matplotlib.pyplot as plt


def GpCalculationStack(gp_im_stack, cell_im_stack, thresh):
    """Generates calculations and plots of the GP values of an organelle within a cell. 

    Parameters
    ----------
    gp_im_stack : file path
        Filepath of the standard image stack.

    cell_im_stack : file path
        Filepath of the cell image stack. 

    thresh : float/int
        Value to set the threshold pixel value when generating the mask 
        organelle. 

    Returns
    -------
    gp_meas : float
        Value of the calculated measure gp.

    gp_corrected : float
        Value of the calculated corrected gp.

    gp_final : float
        Value of the mean final gp.

    final.tif : file
        File saved to device containing the plots of the gp value.
    """

    # finds GP measured value
    im = io.imread(gp_im_stack, plugin="tifffile")  # reads standard images stack
    # im = im.astype(np.int8)

    array_440 = im[:, :, 0]  # indexes the 440 image
    im_440 = np.mean(array_440)  # mean value of 440 image

    array_490 = im[:, :, 1]  # indexes the 490 image
    im_490 = np.mean(array_490)  # mean value of the 490 image

    im_add = np.add(im_440, im_490)  # adds the mean values of the images
    im_sub = np.subtract(im_440, im_490)  # subtracts the mean values of the images
    gp_meas = im_sub / im_add  # divides difference from the sum

    # calculates the GP corrected value
    gp_ref = 0.207  # sets reference value
    gp_numerator = gp_ref + (gp_ref * gp_meas) - gp_meas - 1  # numerator of the corrected equation
    gp_denominator = gp_meas + (gp_ref * gp_meas) - gp_ref - 1  # denominator of the correted equation
    gp_corrected = gp_numerator / gp_denominator  # divides numerator from denominator to give corrected value

    fig = plt.figure(figsize=(30, 30))

    cell_im = io.imread(cell_im_stack, plugin="tifffile")  # reads image stack of the cell

    cell_mask = cell_im[:, :, 2]  # indexes third image as organelle to convert to mask

    # creates mask of 0 and 1 to cover the region of interest (ROI)
    with np.nditer(cell_mask, op_flags=['readwrite']) as it:
        for val in it:
            if val > thresh:
                val[...] = val / val
            else:
                val[...] = 0
    # plots mask
    fig.add_subplot(4, 3, (1))
    plt.imshow(cell_mask, cmap='gray')

    # Calculates GP final with the 440 and 490 equation
    cell_im_440 = cell_im[:, :, 0]  # indexes first image as 440 channel
    cell_im_490 = cell_im[:, :, 1]  # indexes second image as 490 channel

    # plots 440                
    fig.add_subplot(4, 3, (2))
    plt.imshow(cell_im_440, cmap='gray')

    # plots 490                
    fig.add_subplot(4, 3, (3))
    plt.imshow(cell_im_490, cmap='gray')
    gim_490 = (gp_corrected * cell_im_490)

    cell_im_add = cell_im_440 + gim_490  # adds the channels (denominator)

    # plots the addition                  
    fig.add_subplot(4, 3, (4))
    plt.imshow(cell_im_add, cmap='gray')

    cell_im_sub = cell_im_440 - gim_490  # subtracts the channels (numerator)

    # plots the subtraction                 
    fig.add_subplot(4, 3, (5))
    plt.imshow(cell_im_sub, cmap='gray')

    # Heatmap code:
    # sets the divide by zero errors to 200 as placeholders
    cell_im_div = np.divide(cell_im_sub, cell_im_add, out=np.full_like(cell_im_sub, 200),
                            where=cell_im_add != 0, casting='unsafe')

    # sets all remaining zeros to 200 as placeholders
    with np.nditer(cell_im_div, op_flags=['readwrite']) as it:
        for val in it:
            if val == 0:
                val[...] = 200

    # applies the mask to the final image division
    gp_im = np.multiply(cell_im_div, cell_mask)

    # to remove the background, set all 0s to 100
    # remove placeholders and returns values back to 0
    with np.nditer(gp_im, op_flags=['readwrite']) as it:
        for val in it:
            if val == 0:
                val[...] = None
            elif val == 200:
                val[...] = 0

    # plots the division                 
    fig.add_subplot(4, 3, (6))
    plt.imshow(cell_im_div, cmap='gray', vmin=-1, vmax=1)

    # turns array into a data frame
    df = pd.DataFrame(gp_im)

    # creates a heatmap
    fig.add_subplot(4, 3, (7))  # opens a new figure and sets the size
    g = sns.heatmap(df, vmin=-1, vmax=1, cmap='gist_rainbow', xticklabels=False, yticklabels=False)
    g.set_facecolor('xkcd:black')  # sets the background color to black
    g.set_aspect('equal')

    # Histogram code:

    # divides the difference from the sum and sets any division by 0 error to 100 as a placeholder
    cell_gp = np.divide(cell_im_sub, cell_im_add, out=np.full_like(cell_im_sub, 100), where=cell_im_add != 0,
                        casting='unsafe')

    gp_array = np.multiply(cell_gp, cell_mask)  # applies the mask
    gp_condensed = gp_array[gp_array != 0]  # removes any values not in the ROI

    # replaces all placeholders with 0
    with np.nditer(gp_condensed, flags=['zerosize_ok'], op_flags=['readwrite']) as it:
        for val in it:
            if val == 100:
                val[...] = 0
            else:
                val[...] = val

    gp_final = np.mean(gp_condensed)  # final GP value with only values in ROI

    # plots a normalized histogram of these values (area under the histogram = 1)
    fig.add_subplot(4, 3, (8))
    plt.hist(gp_condensed, bins=200, range=[-1, 1], align='mid', alpha=0.5)  # 40 for 0.5 step
    plt.xlabel('GP Value')
    plt.ylabel('Number of pixels')

    df_condensed = pd.DataFrame(gp_condensed)  # creates a data frame of values in ROI
    # plots a PDF curve
    fig.add_subplot(4, 3, (9))
    norm_curve = sns.kdeplot(data=df_condensed, legend=False)
    norm_curve.set_xlim(-1, 1)
    plt.xlabel('GP Value')

    # saves the figure with whichever plots desired (comment out others)
    fig.savefig("final.tif")

    # writes .txt file with GP values
    gp_m = str(gp_meas)
    gp_c = str(gp_corrected)
    gp_f = str(gp_final)
    file = open("GPs.txt", "w")
    file.write(gp_m + ', ' + gp_c + ', ' + gp_f)
    file.close()

    return (gp_meas, gp_corrected, gp_final)


# sample calling code
GpCalculationStack('/Users/riyavarma/Desktop/standard.tif', '/Users/riyavarma/Desktop/cell.tif', 40)

