## ----knit-setup, echo=FALSE, results='hide', eval=TRUE, cache = FALSE----
rm( list = ls() )
library(knitr)
options(fsl.path = "/usr/local/fsl/")
options(fsl.outputtype = "NIFTI_GZ")
library(knitcitations)
library(RefManageR)
cleanbib()
options(citation_format = "pandoc")
bib <- ReadBib("index.bib", check = TRUE)
opts_chunk$set(cache = TRUE, comment = "")
hook1 <- function(x){ gsub("```\n*```r*\n*", "", x) }
hook2 <- function(x){ gsub("```\n+```\n", "", x) }
knit_hooks$set(document = hook2)

## ---- eval=FALSE---------------------------------------------------------
## install.packages("oro.nifti")
## install.packages("fslr")
## install.packages("devtools")
## devtools::install_github("stnava/ITKR") # takes long
## devtools::install_github("stnava/ANTsR") # takes long
## devtools::install_github("muschellij2/extrantsr")

## ---- echo=FALSE---------------------------------------------------------
files = c(t1 = "T1.nii.gz",
          t2 = "T2.nii.gz", 
          flair = "FLAIR.nii.gz",
          pd = "PD.nii.gz")

## ------------------------------------------------------------------------
files

## ----read_t1, message=FALSE----------------------------------------------
library(fslr)
base_t1 = readnii(files["t1"])

## ----ortho, message=FALSE------------------------------------------------
library(fslr)
fslr::ortho2(base_t1)

## ----ortho_overlay_noalpha, message=FALSE--------------------------------
over_50 = mask_img(base_t1, base_t1 > 40)
fslr::ortho2(base_t1, over_50)

## ----image, message=FALSE------------------------------------------------
image(base_t1, z = 55, plot.type = "single")

## ----overlay, message=FALSE----------------------------------------------
over_50[over_50 <= 0] = NA
over_50 = cal_img(over_50)
overlay(base_t1, over_50, z = 55, plot.type = "single")

## ----fsl_bc--------------------------------------------------------------
bc_t1 = fsl_biascorrect(file = base_t1)

## ----bias_corr-----------------------------------------------------------
library(extrantsr)
n4_t1 = bias_correct(file = base_t1, correction = "N4", retimg = TRUE)

## ----bias----------------------------------------------------------------
library(scales)
ratio = finite_img(n4_t1 / base_t1)
ortho2(n4_t1, ratio, col.y = alpha(hotmetal(), 0.5))

## ----bet-----------------------------------------------------------------
ss_t1 = fslbet(n4_t1, outfile = "SS_Image")
ortho2(ss_t1)

## ----ortho_overlay, message=FALSE----------------------------------------
mask = ss_t1 > 0
ortho2(base_t1, y = mask, col.y = alpha("red", 0.5))

## ----cropped_ortho_overlay, message=FALSE--------------------------------
cropped = dropEmptyImageDimensions(ss_t1)
image(cropped, z = floor(dim(cropped)[3]/2), plot.type = "single")

## ----spinning_brain, message=FALSE---------------------------------------
devtools::source_gist("bd40d10afabc503d71e8")

## ----reg_flair-----------------------------------------------------------
reg_flair = flirt(infile = files['flair'], reffile = n4_t1, dof = 6)

## ----ants_reg------------------------------------------------------------
ants_reg_flair = registration(
  filename = files['flair'], 
  template.file = n4_t1, 
  typeofTransform = "Rigid")

## ----double_ortho--------------------------------------------------------
double_ortho(reg_flair, n4_t1)

## ----preproc_within, eval=FALSE------------------------------------------
## proc_images = preprocess_mri_within(
##   files = files[c("t1", "t2", "pd", "flair")],
##   maskfile = ss_t1 > 0)

## ----preproc_across, eval=FALSE------------------------------------------
## outfiles = gsub("[.]nii", '_process.nii', files)
## preprocess_mri_across(
##   baseline_files = files[c("base_t1", "base_t2", "base_pd", "base_flair")],
##   followup_files = files[c("f_t1", "f_t2", "f_pd", "f_flair")],
##   baseline_outfiles = outfiles[c("base_t1", "base_t2", "base_pd", "base_flair")],
##   followup_outfiles = outfiles[c("f_t1", "f_t2", "f_pd", "f_flair")],
##   maskfile = "Brain_Mask.nii.gz")

## ----fast----------------------------------------------------------------
fast_t1 = fast(ss_t1, opts = "--nobias")

## ----otropos-------------------------------------------------------------
tissue_seg = otropos(
  a = ss_t1,
  x = mask)

## ----double_ortho_tissue-------------------------------------------------
double_ortho(ss_t1, tissue_seg$segmentation)

## ----malf----------------------------------------------------------------
# double_ortho(ss_t1, tissue_seg)

