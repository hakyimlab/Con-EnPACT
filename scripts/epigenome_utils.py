import h5py

def query_epigenome(chr_num, center_bp, enfref_dir , num_bins=896, tracks=-1):
    """
    Parameters:
        path_to_enfref (str): path to the directory containing the concatenated reference epigenome files
        chr_num (int/string): chromosome number
        center_bp (int): center base pair position (1-indexed)
        num_bins (int): number of bins to extract centered around center_bp (default: 896) 
            note: if the number of bins is even, the center bin will be in the second half of the array
        tracks (int list): list of tracks to extract (default: all 5313 tracks)

    Returns:
        epigen (np.array): enformer predictions centered at center_bp of shape (num_bins, len(tracks))
    """

    # from position choose center bin
    center_ind = center_bp - 1
    center_bin = center_ind // 128
    
    half_bins = num_bins // 2
    start_bin = center_bin - half_bins
    end_bin = center_bin + half_bins
    if num_bins % 2 != 0: # if num_bins is odd
        end_bin += 1

    with h5py.File(f"{enfref_dir}/chr{chr_num}_cat.h5", "r") as f:
        # get tracks if list provided
        if tracks == -1:
            epigen = f[f'chr{chr_num}'][start_bin:end_bin, :] 
        else:
            epigen = f[f'chr{chr_num}'][start_bin:end_bin, tracks] 

    return epigen



def collect_epigenome(ind, regions, path_to_predictions, collected_preds_dir, start_bin, end_bin):

    print("Collecting epigenome for individual: ",ind)

    predictions_dict_haplo1 = {}
    predictions_dict_haplo2 = {}

    for i,region in enumerate(regions):

        pred_file_h1 = os.path.join(path_to_predictions,ind,"haplotype1",region+"_predictions.h5")
        pred_file_h2 = os.path.join(path_to_predictions,ind,"haplotype2",region+"_predictions.h5")

        if not (os.path.exists(pred_file_h1) and os.path.exists(pred_file_h2)):
            continue
        
        try:
            with h5py.File(pred_file_h1,"r") as f:
                pred_h1 = f[region][:]
                pred_h1 = pred_h1[start_bin:end_bin,:]
                predictions_dict_haplo1[region] = list(pred_h1.mean(axis=0))
            with h5py.File(pred_file_h2,"r") as f:
                pred_h2 = f[region][:]
                pred_h2 = pred_h2[start_bin:end_bin,:]
                predictions_dict_haplo2[region] = list(pred_h2.mean(axis=0))
        except:
            print("Error in reading file: ",ind,region)
            continue

    predictions_df = pd.DataFrame(predictions_dict_haplo1).T
    predictions_df.to_csv(os.path.join(collected_preds_dir, ind+"_haplo1.txt"),sep="\t", header=False)

    predictions_df = pd.DataFrame(predictions_dict_haplo2).T
    predictions_df.to_csv(os.path.join(collected_preds_dir, ind+"_haplo2.txt"),sep="\t", header=False)

    print("Finished collecting epigenome for individual: ",ind)