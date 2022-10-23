FROM ubuntu:22.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget libboost-filesystem1.74-dev libboost-timer1.74-dev \
    libboost-system1.74-dev libboost-date-time1.74-dev libtbb2-dev 

# install itk
#RUN cd /tmp/ \
#    && wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.2.1/InsightToolkit-5.2.1.tar.gz \
#    && cd /opt/ && tar xzvf /tmp/InsightToolkit-5.2.1.tar.gz && cd /opt/InsightToolkit-5.2.1 \
#    && mkdir bin && cd bin && cmake .. && make

RUN apt-get update -qq && apt-get install -yq libinsighttoolkit5-dev

RUN mkdir /pr2mask && cd /pr2mask/ \
    && git clone https://github.com/mmiv-center/pr2nii.git . \
    && echo "Change this string to make this rebuild on docker build" \
    && cmake . && make

ENTRYPOINT [ "/pr2mask/pr2mask" ]
