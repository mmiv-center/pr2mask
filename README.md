# Integration of semantic segmentations into clinical workflows

In semantic segmentation we assign a name like 'tumor' to a region of interest in an image or volume. Individual pixel values can be highlighted in this way and the resulting images allow physicians to evaluate such data. There are several ways to present such data inside clinical systems. All of them are stored as DICOM objects in the archive or VNA.

- As a black and white mask volume: Creating such mask volumes (label fields) is discorraged because they require image fusion and overlays with real data for effective evaluation. They are required as training data for artificial neural networks and they are also required if biomarkers need to be derived from the images.
- Fused image series: Original images as color secondary captures that highlight the mask as a tinted overlay on the original data. Fused images allow for an easy (visual) evaluation of segmentation accuracy. Because they use alpha-blending image intensity inside the mask is also visible. Fused image series cannot easily be used for measurements as they are not coded as original and image viewers might have some limitations if such data is displayed.
- Structured reports: As a DICOM object they keep numeric values with units in the picture archive. Instead of textual images burned into secondary captures all numeric values are stored in a structured way. Such reports can be extracted again from PACS and processed automatically. Structured reports might not in themselves contain image data and thus might not appear in the image archive next to true image objects. The archive will be able to save them. Users can export them again to gain access to structured information.
- Report images: A report blends textual information with key images or derived images similar to a fused image series. Stored as secondary captures these images show up in the image archive and can contain views of regions of interest together with overlays of volume measures or other derived biomarkers. Information presented in a report image are not structured and need to be interpreted by a human/vision A.I.

This package creates all of the above objects to allow for:

- training and finetuning of new A.I. models (mask volumes)
- visual assessment of masks overlayed on primary images (fused image series)
- machine readable biomarker repository (structured reports)
- visual presentation of key derived measures such as tumor volumes (report images)

## Design

Integration of A.I. into a clinical workflow requires support for two steps. i) Creation of training data for image A.I. models and ii) translation of image predictions from trained A.I. models into information suitable for visual assessment and automated data processing.

### Creation of training data for image A.I. models

To support the creation of training data we provide a module for creating segmentation mask from presentation state objects called <i>pr2mask</i>. The presentation state DICOM files (modality PS) contains the polygon vertices for outlines created inside the PACS. All polygons for an image are converted into a mask slice (inside 1, outside value 0) in DICOM format.

If mask slices created by pr2mask are not done every single slice and the segmentation leaves gaps the program <i>MorphologicalContourInterpolation</i> can fill in these gaps and create a volumetric segmentation mask. The module works on the mask DICOM files created by pr2mask and returns a new DICOM mask series where missing slices are filled in. With an option 'FineTuning' these masks can be adjusted to add more detail based on image contrast. This reduces manual segmentation time as a) not all slices need to be manually segmented and b) the polygon segmentations do not need to be that detailed.

#### Presentation state objects to mask volumes

On PACS systems polygonal traces can be stored as DICOM Presentation State objects. In order to get such an object back from the PACS it is not sufficient to 'export' the data from the user interface. Export functionality might not produce a Presentation State object DICOM file. Instead use findscu/movescu to get a copy of the native PACS data.

The data expected by this tool should be a series of DICOM images with an associated Presentation State (PS) file. For example:

```bash
...
MR.1.3.6.1.4.1.45037.f20c2660a1d6c755111a70239b9c98e101fc1d6873170
MR.1.3.6.1.4.1.45037.fe61dff73b8cbee31e79bee6283b55ab11eb62a51c929
PSg.1.2.752.24.7.2440279901.62644.0.481844.0.1664872093
```

The pr2mask tool can be used to convert this list of files into a new output directory:

```bash
# Create mask from presentations state objects using fixed (-u) ids.
# Option '--nobiomarker' switches off the computation of texture and
# shape measures for the regions of interest to speed up computation.
./pr2mask data /tmp/bla -u --verbose --nobiomarker
```


```bash
output
├── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333da_1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333.1.json
├── data.csv
├── fused
│   └── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333.3
│       ├── MR.1.3.6.1.4.1.45037.01876ada7ae594a2199c9498df53d6908dfa8732f3521.dcm
│       ├── MR.1.3.6.1.4.1.45037.07d92f5f5dc1895214c81bdee980e338d81bfe2878573.dcm
│       ├── MR.1.3.6.1.4.1.45037.091a6b29babb817c5ffcc7199d0c0de4cf9d52c155360.dcm
│       ├── MR.1.3.6.1.4.1.45037.1087d261923046a54358b849263399043ee2a9c5de1f5.dcm
│       ├── MR.1.3.6.1.4.1.45037.2a4fc2fa1af7b14fcff5a80655a26779ce2d700c90fb8.dcm
│       ├── MR.1.3.6.1.4.1.45037.505d5aeace2bd9488666ba73cabc24ac5e26387c177a7.dcm
│       ├── MR.1.3.6.1.4.1.45037.52e56a60f0e452a01a757120733fb14e3ed132554170a.dcm
│       ├── MR.1.3.6.1.4.1.45037.5dc2a7d8b8ad61b3416f46f72102136e53fa7d982b4cc.dcm
│       ├── MR.1.3.6.1.4.1.45037.63584f691f1fb48e65255b26c110d0083d4c63e43f4e6.dcm
│       ├── MR.1.3.6.1.4.1.45037.7249e6e1843f52c0342c42ebb6d2ac732501840c6427a.dcm
│       ├── MR.1.3.6.1.4.1.45037.7b636f575748db007cdc7b077ab52f318831e949e73f8.dcm
│       ├── MR.1.3.6.1.4.1.45037.8a7bc4d3310c2568746927d745c538dbd7eb9cba5c8ba.dcm
│       ├── MR.1.3.6.1.4.1.45037.b600a26baed88f9fd04606f9feaa2087017efa99f3174.dcm
│       ├── MR.1.3.6.1.4.1.45037.bf01b9b7563194659c0084704fe1bf5ecf67e3272f1ba.dcm
│       ├── MR.1.3.6.1.4.1.45037.c1c2bb24534b21bbc494e45d6023c9fb61bbb9fa7fa3b.dcm
│       ├── MR.1.3.6.1.4.1.45037.ddf7d12657d8a4dc8a93e49de74d84776c33e8119bc4d.dcm
│       ├── MR.1.3.6.1.4.1.45037.f20c2660a1d6c755111a70239b9c98e101fc1d6873170.dcm
│       └── MR.1.3.6.1.4.1.45037.fe61dff73b8cbee31e79bee6283b55ab11eb62a51c929.dcm
├── images
│   └── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333da
│       ├── MR.1.3.6.1.4.1.45037.01876ada7ae594a2199c9498df53d6908dfa8732f3521.dcm
│       ├── MR.1.3.6.1.4.1.45037.07d92f5f5dc1895214c81bdee980e338d81bfe2878573.dcm
│       ├── MR.1.3.6.1.4.1.45037.091a6b29babb817c5ffcc7199d0c0de4cf9d52c155360.dcm
│       ├── MR.1.3.6.1.4.1.45037.1087d261923046a54358b849263399043ee2a9c5de1f5.dcm
│       ├── MR.1.3.6.1.4.1.45037.2a4fc2fa1af7b14fcff5a80655a26779ce2d700c90fb8.dcm
│       ├── MR.1.3.6.1.4.1.45037.505d5aeace2bd9488666ba73cabc24ac5e26387c177a7.dcm
│       ├── MR.1.3.6.1.4.1.45037.52e56a60f0e452a01a757120733fb14e3ed132554170a.dcm
│       ├── MR.1.3.6.1.4.1.45037.5dc2a7d8b8ad61b3416f46f72102136e53fa7d982b4cc.dcm
│       ├── MR.1.3.6.1.4.1.45037.63584f691f1fb48e65255b26c110d0083d4c63e43f4e6.dcm
│       ├── MR.1.3.6.1.4.1.45037.7249e6e1843f52c0342c42ebb6d2ac732501840c6427a.dcm
│       ├── MR.1.3.6.1.4.1.45037.7b636f575748db007cdc7b077ab52f318831e949e73f8.dcm
│       ├── MR.1.3.6.1.4.1.45037.8a7bc4d3310c2568746927d745c538dbd7eb9cba5c8ba.dcm
│       ├── MR.1.3.6.1.4.1.45037.b600a26baed88f9fd04606f9feaa2087017efa99f3174.dcm
│       ├── MR.1.3.6.1.4.1.45037.bf01b9b7563194659c0084704fe1bf5ecf67e3272f1ba.dcm
│       ├── MR.1.3.6.1.4.1.45037.c1c2bb24534b21bbc494e45d6023c9fb61bbb9fa7fa3b.dcm
│       ├── MR.1.3.6.1.4.1.45037.ddf7d12657d8a4dc8a93e49de74d84776c33e8119bc4d.dcm
│       ├── MR.1.3.6.1.4.1.45037.f20c2660a1d6c755111a70239b9c98e101fc1d6873170.dcm
│       └── MR.1.3.6.1.4.1.45037.fe61dff73b8cbee31e79bee6283b55ab11eb62a51c929.dcm
├── labels
│   └── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333.1
│       ├── MR.1.3.6.1.4.1.45037.01876ada7ae594a2199c9498df53d6908dfa8732f3521.dcm
│       ├── MR.1.3.6.1.4.1.45037.07d92f5f5dc1895214c81bdee980e338d81bfe2878573.dcm
│       ├── MR.1.3.6.1.4.1.45037.091a6b29babb817c5ffcc7199d0c0de4cf9d52c155360.dcm
│       ├── MR.1.3.6.1.4.1.45037.1087d261923046a54358b849263399043ee2a9c5de1f5.dcm
│       ├── MR.1.3.6.1.4.1.45037.2a4fc2fa1af7b14fcff5a80655a26779ce2d700c90fb8.dcm
│       ├── MR.1.3.6.1.4.1.45037.505d5aeace2bd9488666ba73cabc24ac5e26387c177a7.dcm
│       ├── MR.1.3.6.1.4.1.45037.52e56a60f0e452a01a757120733fb14e3ed132554170a.dcm
│       ├── MR.1.3.6.1.4.1.45037.5dc2a7d8b8ad61b3416f46f72102136e53fa7d982b4cc.dcm
│       ├── MR.1.3.6.1.4.1.45037.63584f691f1fb48e65255b26c110d0083d4c63e43f4e6.dcm
│       ├── MR.1.3.6.1.4.1.45037.7249e6e1843f52c0342c42ebb6d2ac732501840c6427a.dcm
│       ├── MR.1.3.6.1.4.1.45037.7b636f575748db007cdc7b077ab52f318831e949e73f8.dcm
│       ├── MR.1.3.6.1.4.1.45037.8a7bc4d3310c2568746927d745c538dbd7eb9cba5c8ba.dcm
│       ├── MR.1.3.6.1.4.1.45037.b600a26baed88f9fd04606f9feaa2087017efa99f3174.dcm
│       ├── MR.1.3.6.1.4.1.45037.bf01b9b7563194659c0084704fe1bf5ecf67e3272f1ba.dcm
│       ├── MR.1.3.6.1.4.1.45037.c1c2bb24534b21bbc494e45d6023c9fb61bbb9fa7fa3b.dcm
│       ├── MR.1.3.6.1.4.1.45037.ddf7d12657d8a4dc8a93e49de74d84776c33e8119bc4d.dcm
│       ├── MR.1.3.6.1.4.1.45037.f20c2660a1d6c755111a70239b9c98e101fc1d6873170.dcm
│       └── MR.1.3.6.1.4.1.45037.fe61dff73b8cbee31e79bee6283b55ab11eb62a51c929.dcm
├── redcap
│   └── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333.1
│       └── output.json
└── reports
    └── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333.1.dcm
```

The 'images' folder contains a copy of the input DICOM images. The 'labels' folder contains a matching list of mask DICOM images as a new series. The matching between images folder and labels folder is stored in the comma-separated file data.csv. If pr2mask is called more than once new rows will be appended to data.csv. The JSON file contains the polygon information parsed from the input presentations state files. The folder fused contains a new secondary capture image series where the mask volume is 'fused' with the individual DICOM images. Regions inside the mask are displayed in color. The 'reports' folder contains a single secondary capture image file with the structure information such as volumes and connected regions shape measures. The 'redcap' folder contains for each image series a json that lists all the output measures in REDCap format. Here an example:

```bash
[
    {
        "field_name": "boundingbox",
        "record_id": "test01",
        "redcap_event_name": "EventName:01",
        "redcap_repeat_instance": 0,
        "redcap_repeat_instrument": "pr2mask",
        "value": "ImageRegion (0x7f9828926aa8)\n  Dimension: 3\n  Index: [69, 63, 14]\n  Size: [105, 118, 2]\n"
    },
...
```

If the program is used continuously converting new polygonal outlines into masks a fixed UID scheme can be switched on using the option '--uid-fixed'. In this mode the program will create the same output UIDs (Series and SOPImage) based on the information of the input series. A PACS system in 'overwrite' mode should replace existing image series with the latest versions generated by pr2mask. This allows for an interactive label refinement workflow.

### Tranlations of image masks into DICOM objects

Semantic segmentation procedures may result in individual masks (black and white) that match with an individual input image series. The <i>imageAndMask2Report</i> program accepts a DICOM image series with a matching DICOM mask series (as created by pr2mask).

```bash
# Generate outputs from image and mask DICOM in /tmp/blarg folder
# Option '-u' ensures that new UIDs are derived from input
# UIDs. This ensures that we can submit such objects several times
# to PACS and they will replace each other (overwrite mode).
./imageAndMask2Report data/input data/mask /tmp/blarg -u
```

All processing results are stored in an output folder with the following structure:

```bash
.
├── 17c7b0b12cf_17c7b0b12cfae8a73650.1.json
├── data.csv
├── fused
│   └── 17c7b0b12cfae8a736509d963df1389b7e0432b0daadac1939a277b10b968.3
├── images
│   └── 17c7b0b12cfae8a736509d963df1389b7e0432b0daadac1939a277b10b968566
├── labels
│   └── 17c7b0b12cfae8a736509d963df1389b7e0432b0daadac1939a277b10b968.1
├── redcap
│   └── 17c7b0b12cfae8a736509d963df1389b7e0432b0daadac1939a277b10b968.1
└── reports
    └── 17c7b0b12cfae8a736509d963df1389b7e0432b0daadac1939a277b10b968.1_0.dcm
```

The structured values for biomarker recovery are contained in the json formatted file. In order to convert them for PACS use

```bash
./json2SR /tmp/blarg/17c7b0b12cf_17c7b0b12cfae8a73650.1.json
```

This will create a new DICOM object in the /tmp/blarg/ folder as a structured report.

### Use partial polygonal outlines for 3D segmentation

The morphological contour interpolation algorithm allows us to create volume filling segmentations from a few polygonal regions (e.g. done every two slices). Create the partial segmentional mask volume first. Use the '-o' option to prevent the expensive computation of biomarkers, the '-u' option will prevent random series and image instance UIDs from being created, e.g. repeating the process will create stable UIDs based on the input. Use MorphologicalContourInterpolation second pointing it to the generated mask volume from pr2mask.


```bash
./pr2mask data/pr2mask_test_liver /tmp/bla -u -o
./MorphologicalContourInterpolation /tmp/bla/labels/1.3.6.1.4.1.45037.ffc60b90cf48bc76cd655d454f2bf8ae6aaf8ebde42.1 \
				    /tmp/output \
				    --image-series /tmp/bla/images/1.3.6.1.4.1.45037.ffc60b90cf48bc76cd655d454f2bf8ae6aaf8ebde4262 \
				    -u
```

As a last step we can use the generated mask and the original image series to create a report that includes volume measures using imageAndMask2Report.

```bash
./imageAndMask2Report /tmp/bla/images/1.3.6.1.4.1.45037.ffc60b90cf48bc76cd655d454f2bf8ae6aaf8ebde4262 /tmp/output/ /tmp/blarg -u -r mosaic
```

The output folder /tmp/blarg/ will contain the output in multiple formats

### Fused image series

Fused image series display the map ontop of the image data in color. If you have a binary mask as an output this will create the fused images.

```bash
./imageAndMask2Fused /tmp/bla/images/1.3.6.1.4.1.45037.ffc60b90cf48bc76cd655d454f2bf8ae6aaf8ebde4262 /tmp/output/ /tmp/blarg -u -i "Fused mask"
```

If your algorithm creates a map that can be used to show the models confidence in the results you can use imageAndMask2Fused to create vote-map images. For example you can train 5 AI models with each 4/5th of the available training data. Summing up the individual binary masks (max value is 5) you create a vote-map.

```bash
./imageAndMask2Fused /tmp/bla/images/1.3.6.1.4.1.45037.ffc60b90cf48bc76cd655d454f2bf8ae6aaf8ebde4262 /tmp/output/ /tmp/blarg --votemapmax 5 --votemapagree 0.5 -u -i "Fused mask"
```

This would create a fused image series with a blue overlay where the vote map reached or exceeded a value of 0.5 * 5 = 2.5. All other voxel with a non-zero value would be displayed in a yellow overlay color.

## Build these modules

This project depends on the boost, libzip and freetype libraries. Install these before starting to build a non-docker version of this program. 

The easiest way to create the executable is by following the instructions in the Dockerfile:

```bash
docker build --no-cache -t pr2mask -f Dockerfile .
docker run --rm -it pr2mask data/ output
```

Note: build in debug mode:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```

If you want to build in macos you might need to specify the location of qt5 with

```bash
cmake . -DCMAKE_PREFIX_PATH=/opt/homebrew/Cellar/qt@5/5.15.10
```

If you update your build system and you are using a source code version of IKT you might need to rebuild IKT (5.4) with

```
cd InsightToolkit
mkdir bin; cd bin
cmake -DModule_MorphologicalContourInterpolation:BOOL=ON ../
```

Delete the CMakeFiles folder and makefile plus CMakeCache.txt.

## Debugging

To create a report you can do:

```bash
./imageAndMask2Report data/input data/mask /tmp/blarg -u
./json2SR /tmp/blarg/*.json
```
