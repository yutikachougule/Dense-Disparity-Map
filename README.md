# Epipolar Geometry and RANSAC Algorithm for Fundamental Matrix Estimation

This project demonstrates the process of detecting Harris corners, applying the RANSAC algorithm, and estimating the fundamental matrix between two images. It also includes steps for calculating dense disparity maps using the epipolar geometry.

## Overview

The main tasks in this code include:
1. **Loading Images**: Reading a dataset of images (e.g., from a directory containing masks).
2. **Harris Corner Detection**: Detecting Harris corner features from the images.
3. **RANSAC Algorithm**: Estimating the fundamental matrix between two images using RANSAC.
4. **Dense Disparity Map**: Computing the dense disparity map using epipolar geometry and the estimated fundamental matrix.
5. **Normalized Cross-Correlation**: Using NCC to match features between two images.

## Prerequisites

Before running this code, make sure you have the following installed:
- MATLAB
- Image Processing Toolbox (for functions like `imfilter`, `imgaussfilt`, etc.)

## Code Structure

### 1. **Reading Image Dataset**

This section sets up directories for the images and reads the first image from the dataset. The images are stored in an `imageDatastore` object, and the number of images in the dataset is calculated.

### 2. **Loading Images into Arrays**

All the images are read and stored into two arrays:
- `gray_scale_images`: Contains grayscale versions of the images.
- `color_images`: Stores the original color images.

### 3. **Harris Corner Detection**

This part detects Harris corner features in each image. The coordinates of the detected corners are stored in the `harris_features` array.

### 4. **RANSAC and Fundamental Matrix Estimation**

In this step, the features from two selected images are matched using normalized cross-correlation (NCC). The RANSAC algorithm is then applied to estimate the fundamental matrix between these two images.

### 5. **Dense Disparity Map Calculation**

This section calculates the dense disparity map using the fundamental matrix estimated in the previous step. It applies epipolar geometry to find disparities across corresponding points in the two images.

### 6. **Helper Functions**

Several helper functions are defined to support the main algorithm:

- **`getHarrisFeatures`**: Detects Harris corners in an image.
- **`getImageFeatures`**: Extracts the features for an image.
- **`getNCCArray`**: Calculates the Normalized Cross-Correlation (NCC) between image patches.
- **`getRANSACInliers`**: Applies the RANSAC algorithm to find inliers and estimate the fundamental matrix.
- **`getMatrixA`**: Constructs the matrix A for RANSAC inlier estimation.

## Function Definitions

### `getHarrisFeatures`

This function detects Harris corner features in a grayscale image. It applies Gaussian smoothing and the Prewitt operator for gradient computation, followed by non-maximum suppression to identify corner points.

### `getImageFeatures`

Extracts the features from the array of Harris corners, selecting features corresponding to a specific image.

### `getNCCArray`

This function calculates the normalized cross-correlation (NCC) score between image patches centered around the detected Harris corners.

### `getRANSACInliers`

The RANSAC algorithm is applied here to find the best set of inliers that fit the fundamental matrix.

### `getMatrixA`

Constructs the matrix A, which is used in the RANSAC process to estimate the fundamental matrix.

## Example Usage

To use this code with your own dataset, set the directories for your images and follow the steps outlined in the code to load the images, detect Harris features, match features between images, estimate the fundamental matrix using RANSAC, and calculate the dense disparity map.

```matlab
% Set directories for the dataset
wall_dir = fullfile("wall/");
scene_dir = fullfile("mask/");

% Load images and detect Harris corners
curr_dir_handle = imageDatastore(scene_dir);
numOfImages = numel(curr_dir_handle.Files);

% Call the function for Harris feature extraction
harris_features = getHarrisFeatures(gray_scale_images(:,:,1), size(gray_scale_images(:,:,1)), 1.5, 21000, 0.04, 1);
```

## Conclusion
This project demonstrates a full pipeline for epipolar geometry analysis using Harris corner detection, the RANSAC algorithm for fundamental matrix estimation, and disparity map computation. This is useful in computer vision tasks such as stereo vision, structure from motion, and 3D reconstruction.

## Notes:
The code assumes images are stored in a directory and can be read using imageDatastore.

Make sure to adjust paths accordingly if using a different dataset.


