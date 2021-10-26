
data_LAND_prep = 'C:\Users\Utente\OneDrive - Fondazione Istituto Italiano Tecnologia\data_LAND_RS_prep'
cd(data_LAND_prep)

set_file_list = dir('*.mat')

figure_dir = fullfile(data_LAND_prep, '_FIG_ASR+Hurst')

for i_subj = 1:length(set_file_list)
    file_name = set_file_list(i_subj).name(1:end-20);

    cd(figure_dir)

    if ~contains(file_name, 'error')
        try
            append_TIFF2TIFF_fun(file_name);
        catch error_tmp
            disp(error_tmp)
        end
    end
end

disp('...DONE!!')


%% ########################################################################
    function append_TIFF2TIFF_fun(file_name)

    %
    %

        %% file_name = 'a4mb6S3l22_01_RestingStateOpen'

        data_LAND_prep = 'C:\Users\Utente\OneDrive - Fondazione Istituto Italiano Tecnologia\data_LAND_RS_prep';
        %cd(data_LAND_prep)


        %% FIGURE type:
        fig_format = '.jpg';
        %fig_format = '.tiff'

        fig_type1 = 'preproc_notchdata_scroll';
        fig_type2 = 'ASR_result';

        fig_type3 = 'preproc_badSamplesChannels';
        fig_type4 = 'preproc_psd';

        fig_type5 = 'preproc_IC_brain_probabilities';
        fig_type6 = 'preproc_ica_1_28';


        %%
        fig_resize_dim = [1236, 2500];
        resize_factor = 1.5;

        A = imread([file_name, '_' fig_type1 , fig_format]); %'.jpg']);
        B = imread([file_name, '_' fig_type2 , fig_format]);
        C = imread([file_name, '_' fig_type3 , fig_format]);
        D = imread([file_name, '_' fig_type4 , fig_format]);
        E = imread([file_name, '_' fig_type5 , fig_format]);
        F = imread([file_name, '_' fig_type6 , fig_format]);

        A = imresize(A,[fig_resize_dim(1)/resize_factor fig_resize_dim(2)]);
        B = imresize(B,[fig_resize_dim(1)/resize_factor fig_resize_dim(2)]);

        C = imresize(C,[fig_resize_dim(1)*resize_factor fig_resize_dim(2)]);
        F = imresize(F,[fig_resize_dim(1)*resize_factor fig_resize_dim(2)]);

        D = imresize(D,[fig_resize_dim(1)*resize_factor fig_resize_dim(2)]);
        E = imresize(E,[fig_resize_dim(1)*resize_factor fig_resize_dim(2)]);
        %D = imresize(D,jpg_resize_dim);
        %E = imresize(E,jpg_resize_dim);

        fig_all = [A, B; C, F; D, E];

        %fig_all = [A; B; C; F; E; D];

        %figure; imshow([A,B])
        figure; imshow(fig_all);

        suptitle(strrep(file_name,'_','  '));

        disp(['...saving: ' file_name])
        cd(data_LAND_prep)
        saveas(gcf, [ file_name '_REPORT'], 'tiff');
        close 
    end