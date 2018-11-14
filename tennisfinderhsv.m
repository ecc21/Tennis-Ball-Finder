image1=im2double(imread('NewBallsLarge.jpg')); %%% failed 1 found 2
image2=im2double(imread('OneBallCornerLarge.jpg')); %% pass
image3=im2double(imread('OneBallLarge.jpg')); %%% fail, false positives
image4=im2double(imread('OneBallLetteringVerticalLarge.jpg')); % pass
image5=im2double(imread('OneBallShadowsLarge.jpg')); %%% bare pass
image6=im2double(imread('OneBallVerticalLarge.jpg')); %pass
image7=im2double(imread('ThreeBallsCloseUpTouching.jpg')); %% failed
image8=im2double(imread('ThreeBallsNetLarge.jpg')); % passed
image9=im2double(imread('ThreeBallsShadowLarge.jpg')); %found 1, failed other due to sticking together
image10=im2double(imread('TwoBallsLetteringLarge.jpg')); % passed
image11=im2double(imread('TwoBallsShadowLarge.jpg')); % failed: false positives
image12=im2double(imread('TwoBallsTouchingVerticalLarge.jpg')); % failed: sticking together
image13=im2double(imread('TwoBallsVerticalLarge.jpg')); % passed

hsv_image = rgb2hsv(image1);
hue_image = hsv_image(:,:,1);
sat_image = hsv_image(:,:,2);
val_image = hsv_image(:,:,3);
%checking hue colors
imshow(hsv_image);
hp = impixelinfo;
set(hp,'Position',[5 1 300 20]);

hue_thresh_low = 0.11;  %Set threshold for being dark. Chosen experimentally
hue_thresh_high = 0.3;  %Set threshold for being dark. Chosen experimentally
hue_thresh = find((hue_image > hue_thresh_low) & (hue_image < hue_thresh_high));

sat_thresh_low = 0.3;  %Set threshold for being dark. Chosen experimentally
sat_thresh_high = 0.5;  %Set threshold for being dark. Chosen experimentally
sat_thresh = find((sat_image > sat_thresh_low) & (sat_image < sat_thresh_high));

val_thresh_low = 0.6;  %Set threshold for being dark. Chosen experimentally
val_thresh_high = 1.1;  %Set threshold for being dark. Chosen experimentally
val_thresh = find((val_image > val_thresh_low) & (val_image < val_thresh_high));

tw=zeros(size(hue_image)); 
tw_hue=-ones(size(hue_image));   %Make up a new array full of -1's that's the same
tw_sat=-ones(size(sat_image));   %Make up a new array full of -1's that's the same
tw_val=-ones(size(val_image));   %Make up a new array full of -1's that's the same

tw_hue(hue_thresh) = 1;
tw_sat(sat_thresh) = 1;
tw_val(val_thresh) = 1;

%thresh = find( ((tw_hue == 1)&(tw_sat == 1)) | ((tw_hue == 1)&(tw_val == 1)) | ((tw_sat == 1)&(tw_val == 1)) );
thresh = find((tw_hue == 1) &(tw_sat == 1) &(tw_val == 1));
tw(thresh) = 1;

tw = bwareaopen(tw,200);

se = strel('disk',25);
tw = imclose(tw, se);

tw = imfill(tw, 'holes');

% code from https://www.mathworks.com/matlabcentral/answers/380687-how-to-smooth-rough-edges-along-a-binary-image
% convolve image and rethreshold to smooth edges
%doesnt really work, since rounding noise creates a lot of false positives
%windowSize = 15;
%kernel = ones(windowSize) / windowSize ^ 2;
%blurryImage = conv2(tw, kernel, 'same');
%binaryImage = blurryImage > 0.5; % Rethreshold
%%% end of code snippit
%tw = binaryImage;

figure; imshow(tw);
title('thresholded image');


% code from https://www.mathworks.com/help/images/identifying-round-objects.html
[B,L] = bwboundaries(tw);
% Display the label matrix and draw each boundary
figure; imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

stats = regionprops(L,'Area','Centroid','MajorAxisLength','MinorAxisLength');

threshold = 0.75;

radiis = [];
centroids = [];
% loop over the boundaries
for k = 1:length(B)

  % obtain (X,Y) boundary coordinates corresponding to label 'k'
  boundary = B{k};

  % compute a simple estimate of the object's perimeter
  delta_sq = diff(boundary).^2;    
  perimeter = sum(sqrt(sum(delta_sq,2)));
  
  % obtain the area calculation corresponding to label 'k'
  area = stats(k).Area;
  
  % compute the roundness metric
  metric = 4*pi*area/perimeter^2;
  
  % display the results
  metric_string = sprintf('%2.2f',metric);

  % mark objects above the threshold with a black circle
  if metric > threshold
    centroid = stats(k).Centroid;
    radii = (mean([stats(k).MajorAxisLength stats(k).MinorAxisLength],2))/2;
    centroids = [centroids ;centroid];
    radiis = [radiis ;radii];
    plot(centroid(1),centroid(2),'ko');
    %%% append centriod and radii
  end
  
  text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
       'FontSize',14,'FontWeight','bold');
  
end

title(['Metrics closer to 1 indicate that ',...
       'the object is approximately round']);
   
centroids
radiis