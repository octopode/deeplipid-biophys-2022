import numpy as np
import os
from skimage import io


def StandardIter(standard_dir, standard_finish):
    """Generates GP correction values of multiple standard images and outputs their mean.

            Parameters
            ----------

            standard_dir : directory path
                Path of the directory with standard images as .tif files.

            standard_finish : directory path
                Path of the directory where you want to save your outputs.

            Returns
            -------

             GPvalues.txt
                Text file containing list of GP correction values in order as well
                as their mean value.

            """
    # sets parameters
    directory = standard_finish
    filename = "GPvalues.txt"
    file_path = os.path.join(directory, filename)  # creates a new file path for the saved images
    gp_file = open(file_path, "w")  # opens a text file for GP final values
    gp_list = [] # opens empty list

    # for loops over folder of image stacks
    for file in os.listdir(standard_dir):
        if file.endswith('.tif'):
            im = os.path.join(standard_dir, file)

            std_im = io.imread(im, plugin="tifffile")  # reads image stack of the cell

            array_440 = std_im[:, :, 0]  # indexes the 440 image
            im_440 = np.mean(array_440)  # mean value of 440 image

            array_490 = std_im[:, :, 1]  # indexes the 490 image
            im_490 = np.mean(array_490)  # mean value of the 490 image

            im_add = np.add(im_440, im_490)  # adds the mean values of the images
            im_sub = np.subtract(im_440, im_490)  # subtracts the mean values of the images
            gp_meas = im_sub / im_add  # divides difference from the sum

            # calculates the GP corrected value
            gp_ref = 0.207  # sets reference value
            gp_numerator = gp_ref + (gp_ref * gp_meas) - gp_meas - 1  # numerator of the corrected equation
            gp_denominator = gp_meas + (gp_ref * gp_meas) - gp_ref - 1  # denominator of the correted equation
            gp_corrected = gp_numerator / gp_denominator  # divides numerator from denominator to give corrected value
            gp_c = str(gp_corrected)
            gp_file.write(gp_c + '\n') # writes file with each corrected GP
            gp_list.append(gp_corrected) # adds value to list every iteration

    gp_mean = np.mean(gp_list) # takes mean of list
    gp_m = str(gp_mean)
    gp_file.write('mean: ' + gp_m) # writes file with mean value
    gp_file.close() # closes file

    return (gp_mean)


StandardIter('/Users/riyavarma/Desktop/standard', '/Users/riyavarma/Desktop/standardfinished')

