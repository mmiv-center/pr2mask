FROM ubuntu:24.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget libboost-filesystem1.83-dev libboost-timer1.83-dev \
    libboost-system1.83-dev libboost-date-time1.83-dev libtbb-dev \
    libboost-program-options1.83-dev libfreetype-dev libxml2-dev zlib1g-dev \
    libzip-dev zipcmp zipmerge ziptool

# install itk
RUN cd /tmp/ \
    && wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.4.2/InsightToolkit-5.4.2.tar.gz \
    && cd /opt/ && tar xzvf /tmp/InsightToolkit-5.4.2.tar.gz && cd /opt/InsightToolkit-5.4.2 \
    && sed -i '22i #include <cstdint>' Modules/Filtering/MathematicalMorphology/include/itkMathematicalMorphologyEnums.h \
    && sed -i "s/uint8_t/std::uint8_t/" Modules/Filtering/MathematicalMorphology/include/itkMathematicalMorphologyEnums.h \
    && mkdir bin && cd bin && cmake -DModule_MorphologicalContourInterpolation:BOOL=ON .. && make -j 4 \
    && rm /tmp/InsightToolkit-5.4.2.tar.gz


# install dcmtk
RUN cd / \
    && git clone https://github.com/DCMTK/dcmtk.git \
    && cd dcmtk && cmake . && make -j4

#RUN apt-get update -qq && apt-get install -yq libinsighttoolkit5-dev

RUN mkdir /pr2mask && cd /pr2mask/ \
    && git clone https://github.com/mmiv-center/pr2mask.git . \
    && echo "Change this string to make this rebuild on docker build" \
    && cmake . && make

ENV REPORT_FONT_PATH=/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf
ENV DCMDICTPATH=/usr/share/libdcmtk16/dicom.dic

ENTRYPOINT [ "/pr2mask/pr2mask" ]
