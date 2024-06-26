FROM ubuntu:22.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget libboost-filesystem1.74-dev libboost-timer1.74-dev \
    libboost-system1.74-dev libboost-date-time1.74-dev libtbb2-dev \
    libboost-program-options1.74-dev libfreetype-dev libxml2-dev zlib1g-dev \
    libzip-dev zipcmp zipmerge ziptool

# install itk
RUN cd /tmp/ \
    && wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.3.0/InsightToolkit-5.3.0.tar.gz \
    && cd /opt/ && tar xzvf /tmp/InsightToolkit-5.3.0.tar.gz && cd /opt/InsightToolkit-5.3.0 \
    && mkdir bin && cd bin && cmake -DModule_MorphologicalContourInterpolation:BOOL=ON .. && make -j 4 \
    && rm /tmp/InsightToolkit-5.3.0.tar.gz

# install dcmtk
RUN cd / \
    && git clone https://github.com/DCMTK/dcmtk.git \
    && cd dcmtk && cmake . && make -j4

#RUN apt-get update -qq && apt-get install -yq libinsighttoolkit5-dev

RUN mkdir /pr2mask && cd /pr2mask/ \
    && git clone https://github.com/mmiv-center/pr2nii.git . \
    && echo "Change this string to make this rebuild on docker build" \
    && cmake . && make

ENV REPORT_FONT_PATH=/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf
ENV DCMDICTPATH=/usr/share/libdcmtk16/dicom.dic

ENTRYPOINT [ "/pr2mask/pr2mask" ]
