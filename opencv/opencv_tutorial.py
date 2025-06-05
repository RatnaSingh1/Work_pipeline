# -*- coding: utf-8 -*-
"""
Created on Fri May 16 11:09:27 2025

@author: ratsa
"""

import os 
import cv2

#read image



img = cv2.imread('photo.jpg')
print(img.shape)


#write
cv2.imwrite(('photo_out.jpg'), img)

#visualize

cv2.imshow('image', img)

# Visualize image
cv2.imshow('image', img)
cv2.waitKey(0)
cv2.destroyAllWindows()




print(resized_img.shape) 
cv2.imshow('resized_img', resized_img)
cv2.waitKey(0)
cv2.destroyAllWindows()

#or


def rescaleFrame(img, scale=1.5):
    width = int(img.shape[1] * scale)
    height = int(img.shape[0] * scale)  # Correct height calculation
    dimensions = (width, height)
    return cv2.resize(img, dimensions, interpolation=cv2.INTER_AREA)  # Correct function name

resized_img = rescaleFrame(img)  # Call the function

cv2.imshow('Resized Image', resized_img)
cv2.waitKey(0)
cv2.destroyAllWindows()
    
 #basic operation resize and recroping
   
resized_img = cv2.resize(img, (1640, 1480))
   
print(img.shape)  
print(resized_img.shape)   
cv2.imshow('resized_img', resized_img)
   
    
   # cropping
   
cropped_img = img[320:660, 420:840] 
cv2.imshow('img', img) 
cv2.imshow('cropped_img', cropped_img) 

cv2.waitKey(0)
cv2.destroyAllWindows() 


# HSV channels, H=hue, S=saturation,V=value

import cv2
import matplotlib.pyplot as plt


# Convert to HSV (make sure this is correct!)
img_hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

# Check if the conversion was successful
print(type(img_hsv))  # Should be numpy.ndarray
print(img_hsv.shape)  # Should be (height, width, 3)

# Split the image into H, S, V components
h, s, v = cv2.split(img_hsv)

# Show the channels
plt.figure(figsize=[20,5])
plt.subplot(141); plt.imshow(h, cmap='gray'); plt.title("H Channel");
plt.subplot(142); plt.imshow(s, cmap='gray'); plt.title("S Channel");
plt.subplot(143); plt.imshow(v, cmap='gray'); plt.title("V Channel");
plt.subplot(144); plt.imshow(cv2.cvtColor(img, cv2.COLOR_BGR2RGB)); plt.title("Original");

plt.show()

# modify the indiviual channels
h_new = h+50
img_NZ_merged = cv2.merge((h_new,s,v))
img_NZ_rgb = cv2.cvtColor(img_NZ_merged, cv2.COLOR_HSV2RGB)
# Show the channels
plt.figure(figsize=[20,5])
plt.subplot(141); plt.imshow(h, cmap='gray'); plt.title("H Channel");
plt.subplot(142); plt.imshow(s, cmap='gray'); plt.title("S Channel");
plt.subplot(143); plt.imshow(v, cmap='gray'); plt.title("V Channel");
plt.subplot(144); plt.imshow(img_NZ_rgb); plt.title("Modified");

# Save the figure before showing it
plt.savefig('all_channels.png', bbox_inches='tight', dpi=300)  # Saves the whole subplot figure
plt.show()

cv2.imwrite(('modified.jpg'), h_new)

#color spaces

img_rgb =   cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    
cv2.imshow('img', img) 
cv2.imshow('img_rgb', img_rgb) 

cv2.waitKey(8000)
cv2.destroyAllWindows()  
    
#blurring
k_size = 9
img_blur=   cv2.blur(img, (k_size, k_size))
img_gaussian_blur =   cv2.GaussianBlur(img, (k_size, k_size), 5)
img_median_blur =   cv2.medianBlur(img, (k_size))
     
cv2.imshow('img', img) 
cv2.imshow('img_blur', img_blur) 
cv2.imshow('img_gaussian_blur r', img_gaussian_blur ) 
cv2.imshow('img_median_blur', img_median_blur) 
cv2.waitKey(0)
cv2.destroyAllWindows()    
    
  #grey scale

  
img_gray =   cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
       
cv2.imshow('img', img) 
cv2.imshow('img_gray', img_gray) 

cv2.waitKey(0)
cv2.destroyAllWindows()   
    
    
img_gray =   cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(img_gray,80,255, cv2.THRESH_BINARY)
cv2.imshow('img', img) 
cv2.imshow('img_gray', img_gray) 
cv2.imshow('thresh', thresh) 
cv2.waitKey(0)
cv2.destroyAllWindows()   
        
 #drawing
line = cv2.line(img, (100,150), (300,450), (0, 0, 255), 5)# strating and end point of line
rect= cv2.rectangle(img, (200, 350), (450, 650), (0, 0, 255), 5) # size, color, thickness
circle = cv2.circle(img, (600, 650), 25, (0, 0, 255), 5) # size, values for center, color, thickness #
cv2.imshow('line', line) 
cv2.imshow('rect', rect)  
cv2.imshow('circle', circle) 
cv2.waitKey(0)
cv2.destroyAllWindows()    
    
 
#image alignment for transformations and compute homography
# Imports
import cv2
import numpy as np
import matplotlib.pyplot as plt
 
    
 # Read reference image
refFilename = "opencv_tutorial/08_Image_Fetaures_and_Alignment/form.jpg"
print("Reading reference image : ", refFilename)
im1 = cv2.imread(refFilename, cv2.IMREAD_COLOR)
im1 = cv2.cvtColor(im1, cv2.COLOR_BGR2RGB)
cv2.imshow("Reading reference image : ", im1) 


# Read image to be aligned
imFilename = "opencv_tutorial/08_Image_Fetaures_and_Alignment/scanned-form.jpg"
print("Reading image to align : ", imFilename)
im2 = cv2.imread(imFilename, cv2.IMREAD_COLOR)
im2 = cv2.cvtColor(im2, cv2.COLOR_BGR2RGB)
cv2.imshow("Reading reference image : ", im2)   
 
    
 # Display Images

plt.figure(figsize=[20,10]); 
plt.subplot(121); plt.axis('off'); plt.imshow(im1); plt.title("Original Form")
plt.subplot(122); plt.axis('off'); plt.imshow(im2); plt.title("Scanned Form")

#find keypoints in both images
    
 # Convert images to grayscale
im1_gray = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
im2_gray = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
  

# Detect ORB features and compute descriptors.
MAX_NUM_FEATURES = 500
orb = cv2.ORB_create(MAX_NUM_FEATURES)
keypoints1, descriptors1 = orb.detectAndCompute(im1_gray, None)
keypoints2, descriptors2 = orb.detectAndCompute(im2_gray, None)

# Display 
im1_display = cv2.drawKeypoints(im1, keypoints1, outImage=np.array([]), color=(255, 0, 0), flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
im2_display = cv2.drawKeypoints(im2, keypoints2, outImage=np.array([]), color=(255, 0, 0), flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

plt.figure(figsize=[20,10])
plt.subplot(121); plt.axis('off'); plt.imshow(im1_display); plt.title("Original Form");
plt.subplot(122); plt.axis('off'); plt.imshow(im2_display); plt.title("Scanned Form"); 



#step, math key points in features
    
# Match features.
matcher = cv2.DescriptorMatcher_create(cv2.DESCRIPTOR_MATCHER_BRUTEFORCE_HAMMING)
matches = matcher.match(descriptors1, descriptors2, None)
  
# Sort matches by score
#matches.sort(key=lambda x: x.distance, reverse=False)


# Convert matches to a list before sorting
matches = list(matches)  # This ensures it's a sortable list
matches.sort(key=lambda x: x.distance, reverse=False)  # Sort matches by distance

# Now matches should be sorted correctly!


# Remove not so good matches
#numGoodMatches = int(len(matches) * 0.1)
#matches = matches[:numGoodMatches]    
    
  # Draw top matches
im_matches = cv2.drawMatches(im1, keypoints1, im2, keypoints2, matches, None)

plt.figure(figsize=[40,10])
plt.imshow(im_matches); plt.axis('off'); plt.title("Original Form");  
    
# step find homography

# Extract location of good matches
points1 = np.zeros((len(matches), 2), dtype=np.float32)
points2 = np.zeros((len(matches), 2), dtype=np.float32)

for i, match in enumerate(matches):
  points1[i, :] = keypoints1[match.queryIdx].pt
  points2[i, :] = keypoints2[match.trainIdx].pt
  
# Find homography
h, mask = cv2.findHomography(points2, points1, cv2.RANSAC)
 
# step wrap image   

# Use homography to warp image
height, width, channels = im1.shape
im2_reg = cv2.warpPerspective(im2, h, (width, height))
# Display results 
plt.figure(figsize=[20,10]); 
plt.subplot(121); plt.imshow(im1); plt.axis('off'); plt.title("Original Form");
plt.subplot(122); plt.imshow(im2_reg); plt.axis('off'); plt.title("Scanned Form");
















#face anomizer
# Load face detection model
face_cascade = cv2.CascadeClassifier(cv2.data.haarcascades + 'haarcascade_frontalface_default.xml')

image = cv2.imread('photo.jpg')
 # Convert to grayscale for better detection
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Detect faces
faces = face_cascade.detectMultiScale(gray, scaleFactor=1.3, minNeighbors=5, minSize=(30, 30))

# Apply blur to detected faces
for (x, y, w, h) in faces:
    face_region = image[y:y+h, x:x+w]  # Extract face region
    face_blur = cv2.GaussianBlur(face_region, (99, 99), 30)  # Apply strong blur
    image[y:y+h, x:x+w] = face_blur  # Replace original face with blurred face

# Save the anonymized image
cv2.imwrite('anonymized_image.jpg', image)

# Display the result
cv2.imshow('Face Anonymizer', image)
cv2.waitKey(0)
cv2.destroyAllWindows()
    
    
 
    
    
    #open video

video = cv2.VideoCapture('video.mp4')

    #visualize

    ret = True
    while ret:
        ret, frame = video.read()
        if ret:
            cv2.imshow('frame', frame)
            cv2.waitKey(10)
            cv2.destroyAllWindows()
      video.release()  
      cv2.destroyAllWindows()
      
      #  or

    # Open video file
    video = cv2.VideoCapture('video.mp4')

    # Define output video settings
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Codec for .mp4 format
    output_video = cv2.VideoWriter('output_video.mp4', fourcc, 30, (640, 480))  # Adjust FPS & size as needed

    ret = True
    while ret:
        ret, frame = video.read()
        
        if ret:
            # Resize frame
            frame_resized = cv2.resize(frame, (640, 480))
            
            # Write frame to output video file
            output_video.write(frame_resized)
            
            # Display frame
            cv2.imshow('Frame', frame_resized)
            if cv2.waitKey(40) & 0xFF == ord('q'):  # Press 'q' to exit
                break

    # Release video objects
    video.release()
    output_video.release()
    cv2.destroyAllWindows()
 
    
 #video from camera
import cv2
import sys    

s= 0  # 0 is defalt camera in system
if len(sys.argv) > 1:
    s = sys.argv[1]
    
source =cv2.VideoCapture(s)    

win_name = 'Camera Preview'
cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)

while cv2.waitKey(1) != 27:  #Escape
    has_frame, frame = source.read()
    if not has_frame:
        break
    cv2.destroyWindow(win_name, frame)
    
source.release()
cv2.destroyWindow(win_name)


#or (this is more right)


# Get video source (default to 0 if no argument provided)
s = int(sys.argv[1]) if len(sys.argv) > 1 else 0

source = cv2.VideoCapture(s)

win_name = 'Camera Preview'
cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)

while True:
    has_frame, frame = source.read()
    if not has_frame:
        break

    cv2.imshow(win_name, frame)

    # Wait for user input, escape key (27) to close
    key = cv2.waitKey(1)
    if key == 27:  # Escape key
        print("Closing camera...")
        break

# Properly release camera & close window
source.release()
cv2.destroyAllWindows()




import numpy as np

# Open the camera (default is 0)
source = cv2.VideoCapture(0)

# Check if the camera opened correctly
if not source.isOpened():
    print("Error: Cannot access camera")
    exit()

frames = []  # List to store HSV frames

for i in range(10):
    has_frame, frame = source.read()
    if not has_frame:
        print(f"Error capturing frame {i}")
        break

    # Convert frame to HSV
    img_hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
    frames.append(img_hsv)

    # Show the current frame
    cv2.imshow("HSV Frame", img_hsv)
    
    # Wait for a short time to allow frame capture
    cv2.waitKey(100)

# Release camera and close windows
source.release()
cv2.destroyAllWindows()

# Now 'frames' contains the 10 HSV images for analysis
print(f"Captured {len(frames)} HSV frames for analysis.")


#or all 0 to 9 frame analyyze and saved in folder captured video

import cv2
import os

# Create a folder to save images
save_path = "captured_frames"
os.makedirs(save_path, exist_ok=True)

# Open the camera (default is 0)
source = cv2.VideoCapture(0)

if not source.isOpened():
    print("Error: Cannot access camera")
    exit()

for i in range(10):  # Capture 10 frames
    has_frame, frame = source.read()
    if not has_frame:
        print(f"Error capturing frame {i}")
        break

    # Convert frame to HSV
    img_hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)

    # Save original frame
    cv2.imwrite(f"{save_path}/frame_{i}_original.jpg", frame)

    # Save HSV frame
    cv2.imwrite(f"{save_path}/frame_{i}_hsv.jpg", img_hsv)

    print(f"Saved frame {i}: original & HSV")

source.release()
cv2.destroyAllWindows()

print("All frames captured and saved successfully!")

#or only 10 frame analyzsed and saved with original

import cv2
import matplotlib.pyplot as plt
import os

# Create a folder to save images
save_path = "captured_analysis"
os.makedirs(save_path, exist_ok=True)

# Open the camera (default is 0)
source = cv2.VideoCapture(0)

if not source.isOpened():
    print("Error: Cannot access camera")
    exit()

frame_count = 0
final_frame = None  # Store the 10th frame

while frame_count < 10:
    has_frame, frame = source.read()
    if not has_frame:
        print(f"Error capturing frame {frame_count}")
        break

    frame_count += 1
    
    # Display each frame briefly
    cv2.imshow("Capturing Frames", frame)
    cv2.waitKey(100)  # Wait 100ms before capturing next frame

    # Store the 10th frame for processing
    if frame_count == 10:
        final_frame = frame

source.release()
cv2.destroyAllWindows()

# Perform HSV analysis on the 10th frame
if final_frame is not None:
    img_hsv = cv2.cvtColor(final_frame, cv2.COLOR_BGR2HSV)
    h, s, v = cv2.split(img_hsv)

    # Create subplot for visualization
    plt.figure(figsize=[20, 5])
    plt.subplot(141); plt.imshow(cv2.cvtColor(final_frame, cv2.COLOR_BGR2RGB)); plt.title("Original")
    plt.subplot(142); plt.imshow(h, cmap='gray'); plt.title("H Channel")
    plt.subplot(143); plt.imshow(s, cmap='gray'); plt.title("S Channel")
    plt.subplot(144); plt.imshow(v, cmap='gray'); plt.title("V Channel")

    # Save subplot as an image
    plt.savefig(f"{save_path}/HSV_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()

    print("HSV analysis saved successfully!")

else:
    print("No frame captured for analysis.")
