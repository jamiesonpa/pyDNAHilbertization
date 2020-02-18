import numpy as np
import scipy.misc as smp
import math
from math import sqrt
import PIL
from PIL import Image
from hilbertcurve.hilbertcurve import HilbertCurve


def create_image(i, j):
  image = Image.new("RGB", (i, j), "white")
  return image

def nearest_square(num):
    return round(sqrt(num))**2

def create_binary_list(seq_file):
    #purine = AG = 1
    #pyrimidine = TC = 0
    binary_list = []
    for base in seq_file:
        if base.upper() == "A":
            binary_list.append([1,0,0])
        elif base.upper() == "G":
            binary_list.append([0,1,0])
        elif base.upper() == "T":
            binary_list.append([0,0,1])
        elif base.upper() == "C":
            binary_list.append([.5,.5,0])

    return binary_list



seq_file = open("DNAexample.txt", "r")
seq_len = len(seq_file.read())
nsqare_seqlen = nearest_square(seq_len)
nsqrt = math.sqrt(nsqare_seqlen)
tuplearray = create_binary_list((open("DNAexample.txt", "r")).read())
pixelarray = []
blankpixel = [1,1,1]
pixelcounter = 0
pixelXposition = 0 
pixelYposition = 0 # we start at the top left and write the image like a book
movingwindowparameter = 2 #this is to allow the user to have control over the window upon which the program integrates base-converted RGB values before and after the pixel currently being iterated on


order_check_counter = 1
while math.pow(2, order_check_counter) < nsqrt:
    order_check_counter = order_check_counter + 1

hilbert_dimension_x = math.pow(2, order_check_counter)
hilbert_order = order_check_counter 
hilbert_dimensions = 2

print("The length of the DNA file is: " + str(seq_len))
print("The nearest square is: " + str(nsqare_seqlen) + ", which has a square root of " + str(nsqrt))
print("Hilbert order is: " + str(hilbert_order))

curve = HilbertCurve(hilbert_order, hilbert_dimensions) 
max_hilbert_length = 2**(hilbert_order*hilbert_dimensions)-1
# img = Image.new('RGB', (round(nsqrt), round(nsqrt)), color = 'white')
img = Image.new('RGB', (round(hilbert_dimension_x), round(hilbert_dimension_x)), color = 'white')
tuplearray = tuplearray[0:max_hilbert_length]
for base in tuplearray:
    # we start at the top left and write the image like a book
    if pixelXposition == nsqrt + 1:
        pixelXposition = 1
        pixelYposition = pixelYposition+1

    # for the first few pixels (1/2 the moving window size) we just fill in with blank pixels)
    if pixelcounter < movingwindowparameter:
        pixelarray.append(blankpixel)
        pixelcounter = pixelcounter + 1

    #once we are past the initial part, we can start integrating the data from before and after the pixel we are currently on
    else:
        pixelbluechannelavg = 0
        pixelredchannelavg = 0
        pixelgreenchannelavg = 0
        pixel_redchannel_array = []
        pixel_greenchannel_array = []
        pixel_bluechannel_array = []


        secondary_pixel_counter = -movingwindowparameter
        while secondary_pixel_counter < movingwindowparameter:
            pixel_redchannel_array.append(tuplearray[pixelcounter + secondary_pixel_counter][0])
            secondary_pixel_counter = secondary_pixel_counter + 1
        secondary_pixel_counter = -movingwindowparameter
        while secondary_pixel_counter < movingwindowparameter:
            pixel_greenchannel_array.append(tuplearray[pixelcounter + secondary_pixel_counter][1])
            secondary_pixel_counter = secondary_pixel_counter + 1
        secondary_pixel_counter = -movingwindowparameter
        while secondary_pixel_counter < movingwindowparameter:
            pixel_bluechannel_array.append(tuplearray[pixelcounter + secondary_pixel_counter][2])
            secondary_pixel_counter = secondary_pixel_counter + 1
        
        redchannel_sum = 0
        bluechannel_sum = 0
        greenchannel_sum = 0
        for redpixel_value in pixel_redchannel_array:
            redchannel_sum = redchannel_sum + redpixel_value
        pixelredchannel_avg = redchannel_sum/len(pixel_redchannel_array)
        for greenpixel_value in pixel_greenchannel_array:
            greenchannel_sum = greenchannel_sum + greenpixel_value
        pixelgreenchannel_avg = greenchannel_sum/len(pixel_greenchannel_array)
        for bluepixel_value in pixel_bluechannel_array:
            bluechannel_sum = bluechannel_sum + bluepixel_value
        pixelbluechannel_avg = bluechannel_sum/len(pixel_bluechannel_array)

        if round(pixelredchannel_avg *256) > 256:
            pixelredchannel_avg = 256
        if round(pixelgreenchannel_avg *256) > 256:
            pixelgreenchannel_avg = 256
        if round(pixelbluechannel_avg *256) > 256:
            pixelbluechannel_avg = 256


        coords = curve.coordinates_from_distance(pixelcounter)
        # print(f'coords(h={pixelcounter}) = {coords}')
        img.putpixel((coords[0],coords[1]), (round(pixelredchannel_avg *256), round(pixelgreenchannel_avg*256), round(pixelbluechannel_avg*256)))


    if pixelcounter+movingwindowparameter + 5 < max_hilbert_length:
        pixelcounter = pixelcounter + 1
        pixelXposition = pixelXposition+1
    else:
        break


img.save('DNAimage.png')
    
# # pixels = new.load()
# #print("The length of the DNA file is: " + str(sequencelength))
# #print("The nearest square is: " + str(nsqareseqlen) + ", which has a square root of " + str(nsqrt))

