# OLAS
Source code for the Overlap-Add Stereogram method. See http://www.computationalimaging.org/publications/holographic-near-eye-displays-based-on-overlap-add-stereograms-siggraph-asia-2019/ for details.

To run, download the light fields from https://drive.google.com/file/d/1qmWAVQQRNAbus2koYrFIGzaphnvhye_t/view and place the contents in the `data/` folder. Other CGH algorithms (found in `utils/`) are included for easy comparison. Uncomment the appropriate lines in `make_SLM_pattern.m` to switch between algorithms.

If you use the code provided in this repository or find the paper linked above relevant to your work, please cite the paper linked above. The BibTeX citation for this code/paper is below:

    @article{Padmanaban:2019:OLAS,
        author = {N. Padmanaban and Y. Peng and G. Wetzstein},
        title = {{Holographic Near-Eye Displays Based on Overlap-Add Stereograms}},
        journal = {ACM Trans. Graph. (SIGGRAPH Asia)},
        issue = {38},
        number = {6},
        year = {2019},
    }
