# Background

This Shiny app includes structural variants that were within +1 Mb of 50 neurodegenerative disease associated genes. The structural variants are mapped and filtered by GATK-SV and include thus variants outside the high-quality subset.

The Shiny app displays data tables and prerendered images of structural variants.

# Usage

1. Download the compressed Shiny app folder from this directory to your computer, do not change internal directory order or create new directories while uncompressing
2. Open app.R
3. Make the shinyapp folder into your working directory 
4. Make sure you have the right R packages and versions to run the app (see session info below)
5. Run the code


# Session info of app

> R version 4.1.3 (2022-03-10).  
> Platform: x86_64-apple-darwin17.0 (64-bit).  
> Running under: macOS Catalina 10.15.7.  
>   
> Matrix products: default.  
> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib.  
> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib.  
>   
> locale:  
> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8. 
>   
> attached base packages:  
> [1] stats     graphics  grDevices utils     datasets  methods   base       
>   
> other attached packages:  
> [1] datamods_1.3.4     shinyWidgets_0.7.4 shiny_1.7.3         
> 
> loaded via a namespace (and not attached):  
>  [1] zip_2.2.2         Rcpp_1.0.9        pillar_1.8.1      bslib_0.4.1       compiler_4.1.3    cellranger_1.1.0  later_1.3.0       jquerylib_0.1.4   
>  [9] forcats_0.5.2     tools_4.1.3       digest_0.6.30     memoise_2.0.1     tibble_3.1.8      jsonlite_1.8.3    lifecycle_1.0.3   pkgconfig_2.0.3   
> [17] rlang_1.0.6       openxlsx_4.2.5.1  cli_3.4.1         rstudioapi_0.14   crosstalk_1.2.0   yaml_2.3.6        curl_4.3.3        haven_2.5.1      
> [25] rio_0.5.29        fastmap_1.1.0     htmlwidgets_1.5.4 sass_0.4.2        vctrs_0.5.0       hms_1.1.2         DT_0.26           glue_1.6.2       
> [33] data.table_1.14.4 R6_2.5.1          fansi_1.0.3       readxl_1.4.1      foreign_0.8-82    magrittr_2.0.3    promises_1.2.0.1  ellipsis_0.3.2   
> [41] htmltools_0.5.3   mime_0.12         xtable_1.8-4      httpuv_1.6.6      utf8_1.2.2        stringi_1.7.8     cachem_1.0.6     
