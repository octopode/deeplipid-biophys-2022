import numpy as np
import pandas as pd
import seaborn as sns
import os
from skimage import io
import matplotlib.pyplot as plt


def GivenGpIter(gp_corrected, cell_im_dir, final_dir, thresh):
    """Uses a previously calculated GP corrected value and applies
        it to a mask of a cell organelle to calculate its GP values.
        Outputs the data in the form of a histogram, heatmap, and
        probability density curve.

        Parameters
        ----------
        gp_corrected : float
            Calculated GP corrected value.

        cell_im_dir : directory path
            Path of the directory where you have saved your cell image stacks.

        final_dir : directory path
            Path of the directory where you want to save your outputs.

        thresh : float/int
            Value to set the threshold pixel value when generating the mask
            organelle.

        Returns
        -------

        heatmap.tif : file
            TIFF file saved containing heatmap data.

        histogram.tif : file
            TIFF file saved containing data represented by histogram.

        PDF.tif : file
            TIFF file saved containing data represented by a probability density curve.

         GPvalues.txt
            Text file containing list of GP values in order.

        """
    # sets up parameters
    global gp_final # global variable
    directory = final_dir
    filename = "GPvalues.txt"
    file_path = os.path.join(directory, filename) # creates a new file path for the saved images
    gp_file = open(file_path, "w") # opens a text file for GP final values

    # if desired to reverse color map use the following code:
    # orig_map = plt.cm.get_cmap('gist_rainbow')
    # reversed_map = orig_map.reversed()

    # for loops over folder of image stacks
    for file in os.listdir(cell_im_dir):
        if file.endswith('.tif'):
            im = os.path.join(cell_im_dir, file)

            cell_im = io.imread(im, plugin="tifffile")  # reads image stack of the cell
            cell_mask = cell_im[:, :, 2]  # indexes third image as organelle to convert to mask

            # creates mask of 0 and 1 to cover the region of interest (ROI)
            with np.nditer(cell_mask, op_flags=['readwrite']) as it:
                for val in it:
                    if val > thresh:
                        val[...] = val / val
                    else:
                        val[...] = 0

            # Calculates GP final with the 440 and 490 equation
            cell_im_440 = cell_im[:, :, 0]  # indexes first image as 440 channel
            cell_im_490 = cell_im[:, :, 1]  # indexes second image as 490 channel

            gim_490 = (gp_corrected * cell_im_490)
            cell_im_add = cell_im_440 + gim_490  # adds the channels (denominator)
            cell_im_sub = cell_im_440 - gim_490  # subtracts the channels (numerator)

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

            # created new specific file names for each generated images and assigns them to
            # folders to be saved
            filename, file_extension = os.path.splitext(file)
            file_1 = filename + 'heatmap' + '.tif'
            file_2 = filename + 'histogram' + '.tif'
            file_3 = filename + 'PDF' + '.tif'
            full_path = os.path.join(final_dir, file_1)
            full_2_path = os.path.join(final_dir, file_2)
            full_3_path = os.path.join(final_dir, file_3)

            # turns array into a data frame
            df = pd.DataFrame(gp_im)

            # creates a heatmap
            fig = plt.figure(1, figsize=(30,30)) # opens a new figure and sets the size
            g = sns.heatmap(df, vmin=-1, vmax=1, cmap='gist_rainbow', xticklabels=False, yticklabels=False)
            g.set_facecolor('xkcd:black') # sets the background color to black
            g.set_aspect('equal') # auto scales the heatmap to proper dimensions
            fig.savefig(full_path) # saved figure to previously defined path
            plt.clf() # clears the plot for next iteration

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

            gp_f = str(gp_final) # converts integer to string to be written in text file
            gp_file.write(gp_f + '\n')  # writes file with GP final

            # plots a normalized histogram of these values (area under the histogram = 1)
            fig_2 = plt.figure(2)
            plt.hist(gp_condensed, bins=200, range=[-1, 1], align='mid', alpha=0.5)  # 40 for 0.5 step
            plt.xlabel('GP Value')
            plt.ylabel('Number of pixels')
            fig_2.savefig(full_2_path)
            plt.clf()

            # PDF code:
            df_condensed = pd.DataFrame(gp_condensed)  # creates a data frame of values in ROI

            # plots a PDF curve
            fig_3 = plt.figure(3)
            norm_curve = sns.kdeplot(data=df_condensed, legend=False)
            norm_curve.set_xlim(-1, 1)
            plt.xlabel('GP Value')
            fig_3.savefig(full_3_path)
            plt.clf()

    gp_file.close()  # closes text file

    return(gp_final)


# sample calling of function
GivenGpIter(0.13067996717031, '/Users/riyavarma/Desktop/heatmap', '/Users/riyavarma/Desktop/heatmap1', 35000)

