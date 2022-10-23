# Presentation state objects to mask volumes

On PACS systems polygonal traces can be stored as DICOM Presentation State objects. In order to get such an object back from the PACS it is not sufficient to 'export' the data from the user interface. Export functionality might not produce a Presentation State object DICOM file. Instead use findscu/movescu to get a copy of the native PACS data.

The data expected by this tool should be a series of DICOM images with an associated Presentation State (PS) file. For example:

```bash
...
MR.1.3.6.1.4.1.45037.f20c2660a1d6c755111a70239b9c98e101fc1d6873170
MR.1.3.6.1.4.1.45037.fe61dff73b8cbee31e79bee6283b55ab11eb62a51c929
PSg.1.2.752.24.7.2440279901.62644.0.481844.0.1664872093
```

The pr2mask tool can be used to convert this folder into a new output directory:

```bash
output
├── 1.3.6.1.4.1.45037.6e5abe67770669399fd07972ecc312e90d0b5e18333da_1.3.6.1.4.1.45037.189384771121937551420690572701550935364.json
├── data.csv
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
└── labels
    └── 1.3.6.1.4.1.45037.189384771121937551420690572701550935364
        ├── MR.1.3.6.1.4.1.45037.01876ada7ae594a2199c9498df53d6908dfa8732f3521.dcm
        ├── MR.1.3.6.1.4.1.45037.07d92f5f5dc1895214c81bdee980e338d81bfe2878573.dcm
        ├── MR.1.3.6.1.4.1.45037.091a6b29babb817c5ffcc7199d0c0de4cf9d52c155360.dcm
        ├── MR.1.3.6.1.4.1.45037.1087d261923046a54358b849263399043ee2a9c5de1f5.dcm
        ├── MR.1.3.6.1.4.1.45037.2a4fc2fa1af7b14fcff5a80655a26779ce2d700c90fb8.dcm
        ├── MR.1.3.6.1.4.1.45037.505d5aeace2bd9488666ba73cabc24ac5e26387c177a7.dcm
        ├── MR.1.3.6.1.4.1.45037.52e56a60f0e452a01a757120733fb14e3ed132554170a.dcm
        ├── MR.1.3.6.1.4.1.45037.5dc2a7d8b8ad61b3416f46f72102136e53fa7d982b4cc.dcm
        ├── MR.1.3.6.1.4.1.45037.63584f691f1fb48e65255b26c110d0083d4c63e43f4e6.dcm
        ├── MR.1.3.6.1.4.1.45037.7249e6e1843f52c0342c42ebb6d2ac732501840c6427a.dcm
        ├── MR.1.3.6.1.4.1.45037.7b636f575748db007cdc7b077ab52f318831e949e73f8.dcm
        ├── MR.1.3.6.1.4.1.45037.8a7bc4d3310c2568746927d745c538dbd7eb9cba5c8ba.dcm
        ├── MR.1.3.6.1.4.1.45037.b600a26baed88f9fd04606f9feaa2087017efa99f3174.dcm
        ├── MR.1.3.6.1.4.1.45037.bf01b9b7563194659c0084704fe1bf5ecf67e3272f1ba.dcm
        ├── MR.1.3.6.1.4.1.45037.c1c2bb24534b21bbc494e45d6023c9fb61bbb9fa7fa3b.dcm
        ├── MR.1.3.6.1.4.1.45037.ddf7d12657d8a4dc8a93e49de74d84776c33e8119bc4d.dcm
        ├── MR.1.3.6.1.4.1.45037.f20c2660a1d6c755111a70239b9c98e101fc1d6873170.dcm
        └── MR.1.3.6.1.4.1.45037.fe61dff73b8cbee31e79bee6283b55ab11eb62a51c929.dcm
```

The 'images' folder in output will contain a copy of the input DICOM images. The 'labels' folder will contain a matching list of mask DICOM images as a new series. The matching between images folder and labels folder is stored in the comma-separated file data.csv. If pr2mask is called more than once new rows will be appended to data.csv. The JSON file contains the polygon information parsed from the input presentations state files.

### Build

The easiest way to create the executable is by following the instructions in the Dockerfile:

```bash
docker build -t pr2mask -f Dockerfile .
docker run --rm -it pr2mask data/ output
```

Note: build in debug mode:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```
