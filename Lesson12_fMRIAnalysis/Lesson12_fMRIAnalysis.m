%%%%%%%%%%%%%%%%%%%%%%%
% Analyzing fMRI data %
%%%%%%%%%%%%%%%%%%%%%%%

% Today we are going to look at fMRI data from a short visual experiment:

% The subject saw videos of coherent motion, incoherent motion, and
% biological motion in separate blocks.


clear
path = 'C:\Research\Courses\Matlab Basics\Lesson12_fMRIAnalysis\OM_Visual\';

% The folder above contains 180 dicom images corresponding to the 180 TRs
% of the experiment. Let's load them and take a look:

% Find the names of all the dicom files in the fMRI folder:
filenames = dir([path '*.dcm']);
% Now let's load them:
for i = 1:180
    fmri(:,:,i) = dicomread([path filenames(i).name]);
end

% Run a movie of the fMRI data to see that it looks good:
clim = [200 1100];
figure()
for tr = 1:180
    imagesc(fmri(:,:,tr),clim);
    colormap(gray)
    pause(0.1)
end

% Instead of looking at the entire mosaic, let's look at a single slice:
clim = [200 1100];
figure()
for tr = 1:180
    imagesc(fmri(161:240,161:240,tr),clim);
    colormap(gray)
    pause(0.1)
end


% Now let's focus on three single voxels and plot their time-courses:
voxel1 = double(squeeze(fmri(160+63,160+24,:)));
voxel2 = double(squeeze(fmri(160+64,160+24,:)));
voxel3 = double(squeeze(fmri(160+69,160+32,:)));

figure();
plot(voxel1);
hold on
plot(voxel2,'g');
plot(voxel3,'r');


% Normalizing fMRI data: percent signal change:

voxel1_persig = 100*(voxel1/mean(voxel1)-1);
voxel2_persig = 100*(voxel2/mean(voxel2)-1);
voxel3_persig = 100*(voxel3/mean(voxel3)-1);

figure();
plot(voxel1_persig);
hold on
plot(voxel2_persig,'g');
plot(voxel3_persig,'r');


corr(voxel1,voxel2)
corr(voxel1,voxel3)

%% Reorganize the data into a 4D matrix that will be easy to access:

rows = [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6];
cols = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6];
for tr = 1:180
    for slice = 1:36
        row_start = 1+(rows(slice)-1)*80; % The "field of view" was at 80*80 during the scan.
        row_end = rows(slice)*80;
        col_start = 1+(cols(slice)-1)*80;
        col_end = cols(slice)*80;        
        fmri_4D(:,:,slice,tr) = squeeze(fmri(row_start:row_end,col_start:col_end,tr));        
    end
end
fmri_4D = double(fmri_4D);

% Make sure the transformation went by ok:
figure()
slice = 20;
for tr = 1:180
    imagesc(squeeze(fmri_4D(:,:,slice,tr)),clim);
    colormap(gray)
    pause(0.1)
end

figure()
slice = 40;
for tr = 1:180
    imagesc(flipdim(squeeze(fmri_4D(:,slice,:,tr))',1),clim);
    colormap(gray)
    pause(0.1)
end

%% Experimental design - building a model/predictor

% The structure of the experiment was as follows:
% 1 = coherent motion, 2 = incoherent motion, 3 = biological motion 
trial_order = [1,3,1,2,3,2,1,1,1,3,2,1,3,3,1,2,1,2,3,2,2,3,3,2];

% Each trial contained 8 seconds of stimulation and 6 seconds of rest
% The movie began with 10 seconds of rest and ended with 14 seconds of
% rest. The TR equaled 2sec.

% Let's create a vector with the same length as the fMRI data (180
% timepoints) that contains the point in time each stimulus was presented:

stim = zeros(1,180);
for i = 1:length(trial_order)
    stim_start = 5 + 7*(i-1);
    stim_end = 8 + 7*(i-1);
    
    stim(1,stim_start:stim_end) = 1;
end

% Add a hemodynamic delay of 2 TRs:

stim_delayed = circshift(stim',2)';


% Plot the stimulus structure before and after shift:
figure();
plot(stim);
hold on
plot(stim_delayed,'g');


%% Correlation analysis - identify voxels that are correlated to the stimulus:

slice = 11;
th = 0.3;

% Turn grayscale slice into an RGB image:
rgb_slice = repmat(squeeze(fmri_4D(:,:,slice,1)),[1 1 3]);
% Scale the values to a range of 0-1:
rgb_slice = rgb_slice/max(max(rgb_slice(:,:,1)));

for x = 1:80
    for y = 1:80
        voxel_tc = squeeze(fmri_4D(x,y,slice,:))';
        if mean(voxel_tc) > 100
            norm_tc = 100 * ((voxel_tc/(mean(voxel_tc)) - 1));
            if corr(norm_tc',stim_delayed') > th
                rgb_slice(x,y,1) = 1;
                rgb_slice(x,y,2) = 0;
                rgb_slice(x,y,3) = 0;
            end
        end
    end
end

% Let's plot the slice:
figure();
image(rgb_slice);


