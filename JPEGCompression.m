clear
close all

%% NOTES
% For our code to reliably work, each the image dimensions must be
% multiples of 8 (e.g. the image must be in one of the standard image 
% aspect ratios, such as 640x480, 1080x720, 1920x1080, etc.). The image can
% be RGB (w/ 3 channels) or grayscale (1 channel).

% In addition, our implementation is not as efficient as real JPEG
% compressors:its runtime scales largely with
% the area of the image. Here, we run tests on images of smaller sizes,
% which are included in this directory, but it may still take a few minutes
% for the code to run. 

% 'quant.csv' must be in the same directory as this code for it to work. 
% The quantization matrix is kept separate for easy access and editing.
%% MAIN

%images already in directory: CTScan1.png, CTScan2.png, dog.jpg
imageFileName = "CTScan1.png"; %name of file of image to be compressed, in the same directory as this file
outputFileName = "CTScan1Compressed.jpg";
qualityFactor = 4; %additional quantization factor by which you want to reduce frequencies of image

%calls main compress method to output DCT + quantization of image
%outputs savings in number storage
%output is in 'pairs.txt'
[X,proportionCompression] = compress(imageFileName,'runLengthEncodingOfCurrentImage.txt','quant.csv',qualityFactor); 

%takes as input the output of compress(), 'pairs.txt'
%then runs run-length decoding on input, then inverse DCT of 
% DCT coefficient matrix and outputs the decompressed image matrix
Y = decompress('runLengthEncodingOfCurrentImage.txt','quant.csv',qualityFactor); 

% visualize and save the final compressed image
disp(proportionCompression)
YY = rescale(Y);

figure(1)
originalImage = imread(imageFileName); 
imshow(originalImage)
title("Original Image")

figure(2)
imshow(YY)
title("Compressed Image")

imwrite(YY,outputFileName,'JPEG');
%% MAIN FUNCTIONS (helper methods in following section)

%compress()
%Takes in an image filename, splits the image into its R,G,B components if they are there,
%then iterates over 8x8 blocks in the image to perform the 2D Discrete 
%Cosine Transform on them, followed by quantization using the 8x8 
%quantization matrix stored in the same directory (quant_matrix_filename).
%This produces a sparse matrix of the DCT transform coefficiencts, 
%which is then transcribed to ouput files using run-length encoding.

%Outputs proportion of numbers needed to be stored from original image to
%with run-length encoding. For example, if original raw image is 256x256x3,
%and the number of pairs stored is 9000, that's 18000 numbers stored vs.
%196608, the compression factor would be ~.09,a ~90% reduction 
%in numbers that need to be stored --

function [compressedMatrix,proportionNumbers] = compress(image_matrix_filename,output_file_name,quant_matrix_filename,quality_factor)
    M = imread(image_matrix_filename);
    numRows = size(M,1);
    numCols = size(M,2);
    numChannels = size(M,3);
    compressedMatrix = zeros([numRows,numCols,numChannels]);
    M = double(M);
    for channel = 1:numChannels
        numHorizontalBlocks = round(numCols/8);
        numVerticalBlocks = round(numRows/8);

        for j = 0:numVerticalBlocks-1
           top = j*8+1;
           bottom = top+7;
           for i = 0:numHorizontalBlocks-1
              left = i*8+1;
              right = left+7;
              block = M(top:bottom,left:right,channel);
              DCT_block = DCT(block);
              DCT_block = quantize(DCT_block,quant_matrix_filename,quality_factor);
              DCT_block = DCT_block + 0;
              compressedMatrix(top:bottom,left:right,channel) = DCT_block;
           end    
        end
    end
    [pairs,pairCount] = runLengthEncoding(compressedMatrix);
    proportionNumbers = pairCount*2/(numRows*numCols*numChannels); 
    writematrix(pairs,output_file_name);
end

%decompress()
%Takes in the run-length encoding file of the quantized DCT cosine coefficients
%outputted by "compress," decodes it to reproduce the original sparse matrix, then
%does reverse quantization and inverse 2D DCT to produce a compressed
%version of the original image

function decompressedMatrix = decompress(compressed_image_matrix_filename,quant_matrix_filename,quality_factor)
    M = runLengthDecoding(compressed_image_matrix_filename);
    numRows = size(M,1);
    numCols = size(M,2);
    numChannels = size(M,3);
    decompressedMatrix = zeros([numRows,numCols,numChannels]);
    
    M = double(M);
    for channel = 1:numChannels
        
        numHorizontalBlocks = round(numCols/8);
        numVerticalBlocks = round(numRows/8);

        for j = 0:numVerticalBlocks-1
           top = j*8+1;
           bottom = top+7;
           for i = 0:numHorizontalBlocks-1
              left = i*8+1;
              right = left+7;
              block = M(top:bottom,left:right,channel);
              inverse_DCT_block = reverse_quantize(block,quant_matrix_filename,quality_factor);
              inverse_DCT_block = inverse_DCT(inverse_DCT_block);
%              disp(inverse_DCT_block)
              decompressedMatrix(top:bottom,left:right,channel) = inverse_DCT_block;
           end    
        end
    end
    decompressedMatrix = round(decompressedMatrix);
end

%% Helper methods

%quantize()
%Takes in M, a 8x8 block of DCT cosine coefficients, and conducts element-wise
%division by the quantization matrix stored at
%"quantization_matrix_filename," followed by an additional division by the
%quality factor. The entries in the resulting matrix are rounded to the
%nearest integer closest to zero and returned.
%This produces the sparse matrix necessary for run-length encoding, 
%by reducing the high index, high frequency cosine coefficients to zero. 

function quantized_matrix = quantize(M,quantization_matrix_filename,quality_factor)
    QM = readmatrix(quantization_matrix_filename);
    quantized_matrix = fix(M./QM./quality_factor);
end


%reverse_quantize()
%Takes in M, a 8x8 block of DCT cosine coefficients, and conducts element-wise
%multiplication by the quantization matrix stored at
%"quantization_matrix_filename," followed by an additional multipication by the
%quality factor. The entries in the resulting matrix are rounded to the
%nearest integer closest to zero and returned. This produces the original DCT 
%cosine coefficients from sparse matrix, except those high index, high frequency
%coefficients that were reduced to 0 in quantize() (as multiplication by
%zero is still zero, making the original coefficients are unrecoverable)

function reverse_quantized_matrix = reverse_quantize(M,quantization_matrix_filename,quality_factor)
    QM = readmatrix(quantization_matrix_filename);
    reverse_quantized_matrix = fix(M.*QM.*quality_factor);
end

%one_dim_dct()
%performs 1-dimensional Discrete cosine transform on a 1xN array x[n] :
%X[k] = [from n = 0 to n = N-1 Σ] x[n]*cos((π/N)*(n+1/2)*k)
%outputs the 1xN array of cosine coefficients

function  X = one_dim_DCT(x)
    N = size(x,2);
    X = [];
    for k = 0:N-1
       sum_n = 0;
       for n = 0:N-1
           sum_n = sum_n + x(1,n+1)*cos((pi/N)*(n+.5)*k);
       end  
       X = [X sum_n];
    end
end    

%DCT()
%performs 1-dimensional Discrete cosine transform on first the rows of 
%a NxN array M, then all of the columns, which is mathematically equivalent
%and more computationally efficient than the 2D DCT.
%outputs the NxN array of cosine coefficients, with higher indices in each
%row and column being the coefficients of the higher frequency cosines.

function M = DCT(M)
    for i = 1:size(M,1)
       M(i,:) = one_dim_DCT(M(i,:)); 
    end
    for j = 1:size(M,2)
       M(:,j) = one_dim_DCT(M(:,j)'); 
    end    
end


%one_dim_inverse_dct()
%performs 1-dimensionalal Inverse Discrete Cosine Transform on a 1xN array X[k]:
%which is an array of cosine coefficients generated from a one dimensionalal
%cosine transform. The formula for inverse DCT:
%f(n) = 2/N * (X[0]/2 + [from k = 1 to k = N-1 Σ] X[k]*cos((π/N)*(n+1/2)*k)
%output is the 1xN array of original values that were used to generate the
%cosine coefficients

function x = one_dim_inverseDCT(X)
    N = size(X,2);
    x = [];
    for n = 0:N-1
        sum_k = 0;
        for k = 1:N-1
            sum_k = sum_k +X(k+1)*cos((pi/N)*k*(n+0.5));
        end
        fn = (2/N)*(0.5*X(1)+sum_k);
        x = [x fn];
    end
end

%inverse_DCT()
%performs 2-dimensional inverse Discrete cosine transform on first the rows of 
%a NxN array M, then all of the columns, which is mathematically equivalent
%and more computationally efficient than the 2D inverse DCT.
%outputs the NxN array of original values used to determine the cosine 
%coefficients that make up the input matrix.

function M = inverse_DCT(M)
    for i = 1:size(M,1)
       M(i,:) = one_dim_inverseDCT(M(i,:)); 
    end
    for j = 1:size(M,2)
       M(:,j) = one_dim_inverseDCT(M(:,j)'); 
    end    
end

% runLengthEncoding()
% Takes in a sparse NxN matrix, reshapes it into a 1xN*N matrix, and does
% run-length encoding on the linear array, e.g. [1 0 0 0 0 2 2 2 1 1 1] 
% becomes [(1,1) (0,4) (2,3) (1,3)]. The method returns an array 
% of these pairs, alongside the original dimensions of the input matrix
% to be printed to a file in compress(). 

function [pairs,pairCount] = runLengthEncoding(M)
    pairs = [size(M,1) size(M,2) size(M,3)];
    for channel = 1:size(M,3)
        MM = reshape(M(:,:,channel),1,[]);
        count = 1;
        pairCount = 0;
        previous_number = MM(1);
        for i = 2:size(MM,2)
           if MM(i) ~= previous_number
               pairs = [pairs previous_number count];
               pairCount = pairCount+1;
               count = 1;
               previous_number = MM(i);
           else
               count = count + 1;
           end
        end
        
        pairs = [pairs MM(i) count];
    end
end

% runLengthDecoding()
% This method takes in the name of the file where the run length encoding,
% i.e. the DCT coefficient matrix dimensions and (number,quantity) pairs,
% is stored, and recasts the pairs to a linear array of DCT coefficients,
% then reshapes them to the NxNxK arrays they originally were, for reverse
% quantization and inverse DCT.

function reshapedM = runLengthDecoding(pairs_text_file)
    M = readmatrix(pairs_text_file);
    
    x_dim = M(1);
    y_dim = M(2);
    z_dim = M(3);

    lineararray = [];
    for i = 4:2:size(M,2)
       number = M(i);
       quantity = M(i+1);
       for j = 1:quantity
           lineararray = [lineararray number];
       end    
    end
    
    reshapedM = reshape(lineararray,[x_dim,y_dim,z_dim]);
end