# FMRI ANALYSIS


rm(list=ls())
library(stringr)
library(plyr)
library(dplyr)
library(fslr)
library(extrantsr)
library(ANTsR)
library(spm12r)
id = "KKI2009-19-"
stub = paste0(id, "fMRI")
filename = paste0(stub, ".nii.gz")

##################################
# Getting the TR 
##################################
x = paste0(stub, ".par")
l = readLines(x)
tr = l  %>% 
    str_subset("epetition time") %>% 
    str_replace("(.*):(.*)", "\\2")  %>% 
    str_trim %>% 
    as.numeric
tr = tr/1000
tr = round(tr, 3)

detrend = function(y) {
    n_timepoints = NROW(new_img)
    timepoints = seq(n_timepoints)
    mod = lm(y ~ timepoints)
    residuals(mod)
}
####################################
# Slice Timing Correction
####################################
# orig_img = readnii(filename)
# new_img = orig_img[30,50,3,]

afilename = paste0("a", filename)
if (file.exists(afilename)){ 
    # aimg = readnii(afilename)
} else {
    aimg = fsl_slicetimer(filename, 
        tr = tr, 
        indexing = "up",
        acq_order = "contiguous",
        outfile = afilename
    )
}
# n4filename = paste0("n4", afilename)
# if (file.exists(n4filename)){ 
#     n4_img = readnii(n4filename)
# } else {
#     n4 = img_ts_to_list(aimg)
#     n4_img = llply(n4, 
#         bias_correct,
#         correction = "N4", .progress = "text")
#     n4_img = img_list_to_ts(n4_img)
#     writenii(n4_img, 
#         filename = n4filename)
# }
# n4_img = aimg

# stat_imgs <- stat_img(n4_img, 
#     func = c("median", "mean"))
# mean_img = stat_imgs$mean
# median_img = stat_imgs$median

# ostat_imgs <- stat_img(aimg, 
#     func = c("median", "mean"))
# omean_img = ostat_imgs$mean
# omedian_img = ostat_imgs$median


stub_filename = afilename
boldImage = check_ants(stub_filename)

avg_img <- getAverageOfTimeSeries(boldImage)

#####################
# Full with Half Max twice the vox size
##################
pdim = voxdim(stub_filename)[1]
fwhm = pdim * 2
numberOfTimePoints <- fslval(
    file = stub_filename, 
        keyword ="dim4")
numberOfTimePoints = as.numeric(
    numberOfTimePoints)

#####################
# Motion Calculation
##################
stub = nii.stub(stub_filename)
moco_file = paste0(stub, 
    "_Motion_Params.rda")
moco_fname = paste0(stub, "_moco_img.nii.gz")

if (all(file.exists(c(moco_file, 
    moco_fname)))) { 
    load(moco_file)
    moco_img = antsImageRead(moco_fname[1])
    motionCorrectionResults$moco_img = 
        moco_img
} else {
    motionCorrectionResults <- 
    antsMotionCalculation(boldImage, 
        fixed = avg_img, 
        moreaccurate = 1,
        txtype = "Rigid")
    save(motionCorrectionResults, 
        file = moco_file)
    moco_img <- 
        motionCorrectionResults$moco_img
    antsImageWrite(moco_img, 
        filename = moco_fname)
}

averageImage <- 
    getAverageOfTimeSeries(moco_img)

moco_params <- 
    motionCorrectionResults$moco_params
moco_params = moco_params %>% 
    select(starts_with("MOCO"))
nuisanceVariables <- moco_params
fd = motionCorrectionResults$fd

###################################
# Motion corrected/realigned image
###################################
maskImage <- getMask(averageImage, 
    mean(averageImage), 
    Inf, cleanup = 2)
averageImage[maskImage == 0] <- 0
boldMatrix <- timeseries2matrix(
    moco_img, 
    maskImage)
DVARS <- computeDVARS(boldMatrix)




anatomical = paste0(id, "MPRAGE.nii.gz")
anat_img = check_ants(anatomical)
oimg = readnii(anatomical)

mask_fname = paste0(
    nii.stub(anatomical),
    "_mask.nii.gz")
if (file.exists(mask_fname)) {
    mask = readnii(mask_fname)
} else {
    bet = fslbet_robust(anat_img)
    mask = bet > 0
    writenii(mask, 
        filename = mask_fname)
}
bet = mask_img(oimg, mask)
bet = bias_correct(bet,
    mask = mask,
    correction = "N4")

get_trans = function(outprefix, 
    warp = TRUE) {
    suff = "0GenericAffine.mat"
    if (warp) {
        suff = c("1Warp.nii.gz",
            suff)
    }
    transforms = paste0(outprefix,
        suff)        
}

# seg = otropos(bet,
#     x = mask, 
#     v = 1)


ccor_file = paste0(stub, 
    "_CompCor.rda")
if (all(file.exists(ccor_file))) { 
    load(ccor_file)
} else {
    highvar <- compcor(
        moco_img, 
        maskImage, 
        ncompcor = 6, 
        variance_extreme = 0.975,
        returnhighvarmatinds = TRUE)
    compCorNuisanceVariables <- compcor(
        moco_img, 
        maskImage, 
        ncompcor = 6, 
        variance_extreme = 0.975)
    save(compCorNuisanceVariables, 
        highvar,
        file = ccor_file)
}

globalSignal <- rowMeans(boldMatrix)

nuisanceVariables <- cbind(
    # globalSignal = globalSignal,
    nuisanceVariables, 
    compCorNuisanceVariables)


nuisanceVariables = 
    scale(nuisanceVariables)
nuisance_mod = 
    lm(boldMatrix ~ nuisanceVariables)
rMatrix <- residuals(nuisance_mod)



rMatrix <- frequencyFilterfMRI(
    rMatrix, 
    tr = tr, 
    freqLo = 0.01, 
    freqHi = 0.1, 
    opt = "trig")

DVARSpostCleaning <- computeDVARS(rMatrix)
# globalSignal <- rowMeans(rMatrix)

cleanBoldImage <- matrix2timeseries(
    boldImage, 
    maskImage, 
    rMatrix)
clean_fname = paste0(stub, 
    "_clean_img.nii.gz")
antsImageWrite(cleanBoldImage, 
    filename = clean_fname)
smoothCleanBoldImage = cleanBoldImage * 1

smoothCleanBoldImage = smoothImage(
    cleanBoldImage, 
    sigma = c(rep(fwhm, 3), 0),
    FWHM = TRUE)

smooth_fname = paste0(stub, 
    "_smooth_clean_img.nii.gz")
antsImageWrite(smoothCleanBoldImage, 
    filename = smooth_fname)

cleanMatrix = timeseries2matrix(
    smoothCleanBoldImage, 
    maskImage)




lev_pca = function(mat, 
    center = TRUE,
    scale = TRUE) {
    mat = scale(mat, center = center, 
        scale = scale)
    
    cov = tcrossprod(mat)
    svd = svd(cov, nv = 0)
    eig = svd$d
    ind_q = eig > mean(eig)    
    print(sum(ind_q))
    U = svd$u[, ind_q]  
    hU = tcrossprod(U)
    lev = diag(hU) 
    return(lev)     
}
lev = lev_pca(cleanMatrix)
med_lev = median(lev)
bad = lev > 3 * med_lev
ind = which.max(lev)
# aaimg = subset_4d(aimg, c(ind, ind +1))




mni_filename = paste0(nii.stub(stub_filename),
    "_to_MNI.nii.gz")
if (file.exists(mni_filename) ) {
    fmri_to_mni = readnii(mni_filename)
} else {
    interpolator = "Linear"
    outprefix = "to_anat"
    trans = get_trans(outprefix)
    if (all(file.exists(trans))) {
        reg = list(fwdtransforms = trans,
            interpolator = interpolator)
    } else {
        reg = registration(
            filename = averageImage, 
            template.file = bet, 
            interpolator = interpolator,
            remove.warp = FALSE,    
            outprefix = outprefix,    
            typeofTransform = "SyN")
    }
    # trans = antsTransformRead(reg$fwdtransforms)
    # mat = matrix(
    #     antsTransformGetParameters(trans),
    #     ncol = 4, 
    #     byrow = FALSE)


    mni_brain = mni_fname(brain = TRUE)
    outprefix = "to_MNI"
    trans = get_trans(outprefix)
    if (all(file.exists(trans))) {
        reg_to_mni = list(fwdtransforms = trans,
            interpolator = interpolator)
    } else {
        reg_to_mni = registration(
            filename = bet, 
            template.file = mni_brain, 
            interpolator = interpolator,
            remove.warp = FALSE,
            outprefix = outprefix,
            typeofTransform = "SyN")
    }

    # orig_to_mni = ants_apply_transforms(
    #     fixed = mni_brain,
    #     moving = oimg, 
    #     interpolator = reg_to_mni$interpolator,
    #     transformlist = reg_to_mni$fwdtransforms)

    tlist = c(reg$fwdtransforms,
        reg_to_mni$fwdtransforms)

    abet = oro2ants(bet)
    amni = check_ants(mni_brain)
    ind = 1
    arr = array(NA, 
        dim = c(dim(amni), numberOfTimePoints))
    for (ind in seq(numberOfTimePoints)) {
        run_img = subset_4d(smoothCleanBoldImage,
            ind = ind)
        fmri_to_anat = antsApplyTransforms(
            fixed = abet,
            moving = run_img, 
            interpolator = reg$interpolator,    
            transformlist = reg$fwdtransforms
            )
        fmri_to_mni = antsApplyTransforms(
            fixed = amni,
            moving = fmri_to_anat, 
            interpolator = reg_to_mni$interpolator,    
            transformlist = reg_to_mni$fwdtransforms)
        arr[,,,ind] = as.array(fmri_to_mni)
        print(ind)
    }
    fmri_to_mni = as.antsImage(arr)
    antsCopyImageInfo(
        reference = smoothCleanBoldImage, 
        target = fmri_to_mni)
    antsImageWrite(fmri_to_mni, 
        filename = mni_filename)    
    #         return(fmri_to_mni)
    #         }, .progress = "text")
    # fmri_to_mni = img_list_to_ts(fmri_to_mni)
    # writenii(fmri_to_mni, 
    #     filename = mni_filename)
}

fmri_to_mni = check_nifti(fmri_to_mni)