%% Read Image Dataset


% Set data set directories
wall_dir = fullfile("wall/");
scene_dir = fullfile("mask/");
%check_board_dir = fullfile("CheckerBoard/");

% Select Directory Handle to be used
curr_dir_handle = imageDatastore(scene_dir);
I = readimage(curr_dir_handle,1);
imshow(I);

% Get number of images in the data set
numOfImages = numel(curr_dir_handle.Files)

%% Load Images into Arrays



% Generating Arrays of original image and Grayscale Images
gray_scale_images = zeros(size(I,1),size(I,2),numOfImages);
color_images = zeros(size(I,1),size(I,2),3,numOfImages);

for i=1:numOfImages
    I = readimage(curr_dir_handle,i);
    grayImage = im2gray(I);
    gray_scale_images(:,:,i) = grayImage;
    color_images(:,:,:,i) = I;
end
disp(numOfImages + " Images stored in 3D array : gray_scale_images")

%% Detect Harris Corner Features


% Array for storing detected feature coordinates
harris_features = [];

for i=1:numOfImages
    [x,y] = getHarrisFeatures(gray_scale_images(:,:,i),size(gray_scale_images(:,:,i)), 1.5, 21000, 0.04,1);
    harris_features = [harris_features; [x,y,zeros(numel(x),1) + i]];
end


%% Implement RANSAC algorithm to find Fundamental Matrix



img1 = 1;
img2 = 2;

img1_features = getImageFeatures(harris_features,img1);
img2_features = getImageFeatures(harris_features,img2);

window_size = 11;
    
% Implement Normatlized Cross-corelation

NCC_ARRAY = getNCCArray(gray_scale_images,img1,img2,getImageFeatures(harris_features,img1),getImageFeatures(harris_features,img2), window_size);

% Homography Estimation and RANSAC
img1_matched = img1_features(NCC_ARRAY(:,1),:);
img2_matched = img2_features(NCC_ARRAY(:,2),:);

img1_matched = fliplr(img1_matched);
img2_matched = fliplr(img2_matched);

    
[inliers, img1_matched, img2_matched] = getRANSACInliers(img1_features,img2_features,NCC_ARRAY);

figure;
% Show matched features - ALL
showMatchedFeatures(gray_scale_images(:,:,img1),gray_scale_images(:,:,img2),img1_matched,img2_matched,"montage");

figure;
% Show matched features with inliers
showMatchedFeatures(gray_scale_images(:,:,img1),gray_scale_images(:,:,img2),img1_matched(inliers,:),img2_matched(inliers,:),"montage", PlotOptions={'o', '+', 'g'});

% Calculate matrix A using only inliers
A_inliers = getMatrixA(numel(inliers),img1_matched(inliers,:),img2_matched(inliers,:),0);

% Calculate F 
[U,S,V] = svd(A_inliers' * A_inliers);
F = reshape(U(:,9),[3 3])';

% Enforce rank-2 constraint on F by zeroing out the smallest singular value
[Uf,Sf,Vf] = svd(F);
Sf(end) = 0;
estimated_fundamental_matrix = Uf*Sf*Vf'

%% Dense Disparity Map    

epipolar_line_params = zeros(3);
% 
% d_window_size = 5
% d_padding = floor(d_window_size/2)
% I2_row = zeros(size(gray_scale_images(:,:,1)));
% 
% for i=1+d_padding:size(gray_scale_images,1)-d_padding
%     for j = 1+d_padding:size(gray_scale_images,1)-d_padding
%         epipolar_line_params = [i,j,1] * estimated_fundamental_matrix;
%         I2_row(j,i) = uint32(-epipolar_line_params(3)/epipolar_line_params(2));
% 
%     end
% end

d_window_size = 3;
d_padding = floor(d_window_size/2);
I2_row = zeros(size(gray_scale_images(:,:,1)));

% Create a meshgrid of coordinates
[x, y] = meshgrid(1:size(gray_scale_images,2), 1:size(gray_scale_images,1));
    

% Reshape the coordinates into a 3 x N matrix
coords = [reshape(x, [], 1)'; reshape(y, [], 1)'; ones(1, numel(x))];

% Compute the epipolar lines using matrix multiplication
epipolar_lines = estimated_fundamental_matrix * coords;

I2_row(:) = uint32(-epipolar_lines(3,:) ./ epipolar_lines(2,:));



for i=1+d_padding:size(gray_scale_images,1)-d_padding
    for j = 1+d_padding:size(gray_scale_images,2)-d_padding
        I1 = gray_scale_images(i-d_padding:i+d_padding,j-d_padding:j+d_padding,1);
            if j < size(gray_scale_images(:,:,2)) - 20
                for k = j:j+20
                    I2 = gray_scale_images(I2_row(i,j)-d_padding:I2_row(i,j)+d_padding, k-d_padding:k+d_padding,2);
                    %NCC = corr2(I1,I2);
                    sumabs(I2 - I1);
    
                end
            end
        
    end

end




%% Function Definations

function [x,y] = getHarrisFeatures(image,image_size,sigma,thresh,k,disp)
    
    % filtering the original image 
    image = imgaussfilt(image,sigma);
    
    % Prewitt OPerator
    Ix = [-1 0 1];
    Iy = [-1 0 1]';
    
    % Convolution
    Gx = imfilter(image,Ix,"symmetric");
    Gy = imfilter(image,Iy,"symmetric");
    
    Gx_squared = Gx.^2;
    Gy_squared = Gy.^2;
    Gx_Gy = Gx .* Gy;

    Gx_squared = imgaussfilt(Gx_squared,sigma);
    Gy_squared = imgaussfilt(Gy_squared,sigma);
    Gx_Gy = imgaussfilt(Gx_Gy,sigma);
    
    % Computing the response function R
    R = zeros(image_size);
    for i=1:image_size(1)
        for j=1:image_size(2)
            A_harris = [Gx_squared(i,j) Gx_Gy(i,j);Gx_Gy(i,j) Gy_squared(i,j)];
            R(i,j) = det(A_harris) - (k*(trace(A_harris))^2);
        end
    end
    
    % Thresholding the R function
    R = (R>thresh) .* R;
    
    % Performing Non Max Supression
    R = imregionalmax(R);
    
    % Finding the pixel coordinates of the non-zero values in th R matrix
    [x,y] = find(R);
    
    % Displaying Image with Harris Features
    if disp == 1
        figure;
        imshow(image,[]);
        hold on;
        plot(y,x,'+g',"LineWidth",1);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for extracting image features 

function image_features = getImageFeatures(image, image_number)
    img = image .* (image_number == image(:,3));
    img = nonzeros(img');
    image_features = reshape(img',[3 size(img,1)/3])';
    image_features = image_features(:,1:2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NCC_ARRAY = getNCCArray(gray_scale_images,img_1,img_2,img1_features,img2_features,window_size)
    
    % Define Image 1 and Image 2
    img1 = gray_scale_images(:,:,img_1);
    img2 = gray_scale_images(:,:,img_2);

    % Padding for consttructing neighbourhood
    padding = floor(window_size/2);
    
    % Defining arrays
    NCC_ARRAY = zeros(1,3);
    final_ncc_score = zeros(1,3);

    for i=1:size(img1_features,1)
    
            prev_score = -1;
    
            % NEIGHBOURHODD in IMAGE 1
            try
                I = img1(img1_features(i,1)-padding:img1_features(i,1)+padding,img1_features(i,2)-padding:img1_features(i,2)+padding);
            catch 
                continue
            end
    
            for j=1:size(img2_features,1)
               
                % NEIGHBOURHOOD in IMAGE 2
                try
                    J = img2(img2_features(j,1)-padding:img2_features(j,1)+padding,img2_features(j,2)-padding:img2_features(j,2)+padding);
                catch
                    continue
                end
    
                % Finding the cross-corelation value
                NCC_score =  corr2(I,J); 
                
                % For ensuring one to one mapping
                if NCC_score>prev_score 
                    final_ncc_score = [i j NCC_score];
                    prev_score = NCC_score;
                end
            end
            
            % Thresholding the value of NCC
            if(final_ncc_score(1,3)>0.98)
                NCC_ARRAY(end+1, :) = final_ncc_score;
            end
    
    end
    
    NCC_ARRAY = NCC_ARRAY(2:end,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [highest_inliers, img1_matched, img2_matched] = getRANSACInliers(img1_features,img2_features,NCC_ARRAY)
    img1_matched = img1_features(NCC_ARRAY(:,1),:);
    img2_matched = img2_features(NCC_ARRAY(:,2),:);
    
    img1_matched = fliplr(img1_matched);
    img2_matched = fliplr(img2_matched);
    
    % Sample random points
    
    num_of_points = 8;
    
    inliers = [];
    highest_inliers = [];
    
    for j=1:100
        try
            A = getMatrixA(num_of_points,img1_matched,img2_matched,1);
        catch
            continue
        end


        [U,S,V] = svd(A' * A);
        F = reshape(U(:,9),[3 3])';
        
        
        % Enforce rank-2 constraint on F by zeroing out the smallest singular value
        [Uf,Sf,Vf] = svd(F);
        Sf(end) = 0;
        F = Uf*Sf*Vf';
    
        inliers = [];
    
        for i= 1:size(img1_matched)
           img2_res = [img1_matched(i,1) img1_matched(i,2) 1] * F * [img2_matched(i,1) img2_matched(i,2) 1]';
           

           if abs(img2_res) < 10e-10
                img2_res;
                inliers(end+1) = i; 
           end
           
           if i == 8
               img2_res;
               [img1_matched(i,1) img1_matched(i,2) 1];
               [img2_matched(i,1) img2_matched(i,2) 1];
               img2_res;
           end
        end

        if numel(inliers)>numel(highest_inliers)
            highest_inliers = inliers
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = getMatrixA(num_of_points,img1_matched,img2_matched,random_sample)
    A = zeros(1,9);

    if random_sample
        random_index = randperm(size(img1_matched,1),num_of_points);
    else
        random_index = 1:num_of_points
    end

    for i=1:num_of_points
        P1 = img1_matched(random_index(i),:);
        P2 = img2_matched(random_index(i),:);
        
        x1 = P1(1);
        x2 = P2(1);
        y1 = P1(2);
        y2 = P2(2);
    
        A(end+1,:) = [x1*x2 x1*y2 x1 y1*x2 y1*y2 y1 x2 y2 1];

    end  
    A = A(2:end,:);
end