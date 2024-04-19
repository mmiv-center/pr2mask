#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"

#include "itkMinimumMaximumImageCalculator.h"
#include "itkScalarImageToHistogramGenerator.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"


#include "gdcmAnonymizer.h"
#include "gdcmAttribute.h"
#include "gdcmDataSetHelper.h"
#include "gdcmDirectoryHelper.h"
#include "gdcmFileDerivation.h"
#include "gdcmFileExplicitFilter.h"
#include "gdcmGlobal.h"
#include "gdcmImageApplyLookupTable.h"
#include "gdcmImageChangePlanarConfiguration.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmImageHelper.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmMediaStorage.h"
#include "gdcmReader.h"
#include "gdcmRescaler.h"
#include "gdcmStringFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkConstantPadImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkScalarImageToRunLengthMatrixFilter.h"

#include "itkContinuousIndex.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkGDCMImageIO.h"

#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/interpolators/catmull_rom.hpp>
#include <codecvt>
#include <locale> // wstring_convert
#include <map>
#include <zip.h>

#include "mytypes.h"
#include "report.h"

using json = nlohmann::json;
using namespace boost::filesystem;

json resultJSON;
bool verbose = false;

using ImageType2D = itk::Image<PixelType, 2>;
using MaskImageType2D = itk::Image<PixelType, 2>;

using CPixelType = itk::RGBPixel<unsigned char>;
using CImageType = itk::Image<CPixelType, 2>;
using LabelType = unsigned short;
using ShapeLabelObjectType = itk::ShapeLabelObject<LabelType, 3>;
using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;
typedef itk::Image<PixelType, 3> ImageType3D;

struct generateImageReturn
{
     CImageType::Pointer keyImage;
     std::vector< std::array<int, 2> > pos;
     std::vector< std::string > text;
     std::vector<int> roi_order;
};

// generate a key image based on an input image mask and the ground truth image
// This key image generator shall generate a mosaic of images, one for each lesion.
generateImageReturn generateKeyImageMosaic(ImageType3D::Pointer image, LabelMapType *labelMap, std::vector<int> resolution, float lowerT, float upperT) {
  if (verbose) {
    fprintf(stdout, "Start generating a key image...\n");
  }
  std::vector<std::vector<float>> labelColors2 = {{0, 0, 0}, {166,206,227}, {31,120,180}, {178,223,138}, {51,160,44}, {251,154,153}, {227,26,28}, {253,191,111}, {255,127,0}, {202,178,214}, {106,61,154}, {255,255,153}, {177,89,40}};

  generateImageReturn returns; // store keyImage and the location and text that should be presented ontop of the image

  // based on the number of labels we will create an image mosaic (MPR for each region of interest)
  int numObjects = labelMap->GetNumberOfLabelObjects();
  if (numObjects < 1) {
    fprintf(stderr, "Error: no object found, refuse to create a key image\n");
    return(returns);
  }

  // lets use the base image resolution of
  int base_image_sizeLW = 512; // we will put an axial square image next to coronal and sagittal on top of each other
  int base_image_sizeLH = base_image_sizeLW;
  int base_image_sizeRW = floor(base_image_sizeLW / 1.618);
  int base_image_sizeRH = floor(base_image_sizeLH / 2.0);

  // the total size of the image is now
  resolution[0] = base_image_sizeLW + base_image_sizeRW;
  resolution[1] = numObjects * base_image_sizeLH;

  // create a new RGB image
  CImageType::Pointer keyImage = CImageType::New();
  returns.keyImage = keyImage;

  using RegionType = itk::ImageRegion<2>;
  RegionType::SizeType size;
  size[0] = resolution[0]; // 512
  size[1] = resolution[1]; // 512

  RegionType::IndexType index;
  index.Fill(0);

  RegionType region(index, size);
  keyImage->SetRegions(region);

  keyImage->Allocate();
  keyImage->FillBuffer(itk::NumericTraits<CPixelType>::Zero);

  RegionType outputRegion = keyImage->GetLargestPossibleRegion();
  //itk::ImageRegionIteratorWithIndex<CImageType> outputRGBIterator(keyImage, outputRegion);

  // 
  // Compute the independently blurred channels for the fused image and all colors in 3D
  // 

  //
  // compute optimal window level for whole image, we will use the values for the curved slice
  //
  using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<ImageType3D>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();
  int minGray = imageCalculatorFilter->GetMinimum();
  int maxGray = imageCalculatorFilter->GetMaximum();

  using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImageType3D>;
  HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(image);
  int histogramSize = 1024;
  histogramGenerator->SetNumberOfBins(histogramSize);
  histogramGenerator->SetHistogramMin(minGray);
  histogramGenerator->SetHistogramMax(maxGray);
  histogramGenerator->Compute();
  using HistogramType = HistogramGeneratorType::HistogramType;
  const HistogramType *histogram = histogramGenerator->GetOutput();
  // lowerT and upperT in percentages 0..1 are set in calling function
  double t1 = -1;
  double t2 = -1;
  double sum = 0;
  double total = 0;
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    total += histogram->GetFrequency(bin, 0);
  }
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    double f = histogram->GetFrequency(bin, 0) / total;
    // fprintf(stdout, "bin %d, value is %f\n", bin, f);
    sum += f;
    if (t1 == -1 && sum > lowerT) {
      t1 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
    }
    if (t2 == -1 && sum > upperT) {
      t2 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
      break;
    }
  }
  //if (verbose) {
  //  fprintf(stdout, "calculated best threshold low: %f, high: %f\n", t1, t2);
  //}

  // we need three 3D images for the red green and blue channel
  // so we can blurr them before using them in the output image
  typedef float FPixelType;
  // typedef itk::Image<FPixelType, 2> FloatImageType;
  using FloatImageType = itk::Image<FPixelType, 3>;
  FloatImageType::Pointer red_channel = FloatImageType::New();
  FloatImageType::Pointer green_channel = FloatImageType::New();
  FloatImageType::Pointer blue_channel = FloatImageType::New();
  ImageType3D::RegionType fusedRegion = image->GetLargestPossibleRegion();
  red_channel->SetRegions(fusedRegion);
  red_channel->Allocate();
  red_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  red_channel->SetOrigin(image->GetOrigin());
  red_channel->SetSpacing(image->GetSpacing());
  red_channel->SetDirection(image->GetDirection());
  green_channel->SetRegions(fusedRegion);
  green_channel->Allocate();
  green_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  green_channel->SetOrigin(image->GetOrigin());
  green_channel->SetSpacing(image->GetSpacing());
  green_channel->SetDirection(image->GetDirection());
  blue_channel->SetRegions(fusedRegion);
  blue_channel->Allocate();
  blue_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  blue_channel->SetOrigin(image->GetOrigin());
  blue_channel->SetSpacing(image->GetSpacing());
  blue_channel->SetDirection(image->GetDirection());

  std::vector<std::pair<int,int>> sizeByROI;
  for (unsigned int roi = 0; roi < labelMap->GetNumberOfLabelObjects(); ++roi) {
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(roi); // the label number is the connected component number - not the one label as mask
    sizeByROI.push_back(std::make_pair(labelObject->GetNumberOfPixels(), roi));
  }
  std::sort(sizeByROI.begin(), sizeByROI.end(), [](auto &left, auto &right) {
    return left.first > right.first; // largest roi first
  });
  // add info as roi_order for later
  //for (int roi_idx = 0; roi_idx < labelMap->GetNumberOfLabelObjects(); ++roi_idx) {
  //  returns.roi_order.push_back(sizeByROI[roi_idx].second);
  //}

  // here we need to compute using the label order, later when we look at the results we will use the returned roi_order?
  // but we don't need both!!
  for (unsigned int n_idx = 0; n_idx < labelMap->GetNumberOfLabelObjects(); n_idx++) {
    unsigned int n = sizeByROI[n_idx].second;
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n); // the label number is the connected component number - not the one label as mask
    int label = labelObject->GetLabel(); // it might be that all regions have the same label and only n is a good choice for the color

    // color is
    std::vector<float> col = labelColors2[ (label % (labelColors2.size()-1)) +1 ];
    itk::Index<3U> index;
    for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
      index = labelObject->GetIndex(pixelId);
      // set this position in all three color images, values are between 0 and 1
      red_channel->SetPixel(index, col[0]/255.0); 
      green_channel->SetPixel(index, col[1]/255.0);
      blue_channel->SetPixel(index, col[2]/255.0);
    }
  }

  // smooth all three channels independently from each other
  using GFilterType = itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>;
  auto gaussFilterR = GFilterType::New();
  gaussFilterR->SetInput(red_channel);
  gaussFilterR->SetVariance(1.5f);
  gaussFilterR->Update();
  FloatImageType::Pointer smoothRed = gaussFilterR->GetOutput();

  auto gaussFilterG = GFilterType::New();
  gaussFilterG->SetInput(green_channel);
  gaussFilterG->SetVariance(1.5f);
  gaussFilterG->Update();
  FloatImageType::Pointer smoothGreen = gaussFilterG->GetOutput();

  auto gaussFilterB = GFilterType::New();
  gaussFilterB->SetInput(blue_channel);
  gaussFilterB->SetVariance(1.5f);
  gaussFilterB->Update();
  FloatImageType::Pointer smoothBlue = gaussFilterB->GetOutput();

  itk::ImageRegionIterator<FloatImageType> redSIterator(smoothRed, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> greenSIterator(smoothGreen, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> blueSIterator(smoothBlue, fusedRegion);

  itk::LinearInterpolateImageFunction<ImageType3D, double>::Pointer interpolator =
    itk::LinearInterpolateImageFunction<ImageType3D, double>::New();
  interpolator->SetInputImage(image);

  itk::LinearInterpolateImageFunction<FloatImageType, double>::Pointer interpolatorRed =
    itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolatorRed->SetInputImage(smoothRed);

  itk::LinearInterpolateImageFunction<FloatImageType, double>::Pointer interpolatorGreen =
    itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolatorGreen->SetInputImage(smoothGreen);

  itk::LinearInterpolateImageFunction<FloatImageType, double>::Pointer interpolatorBlue =
    itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolatorBlue->SetInputImage(smoothBlue);

  float f = 0.4; // weight of the underlay, 0.1 is mostly mask visible

  // TODO: sort roi's by size and start with the largest region of interest (top of the image)
  for (unsigned int roi_idx = 0; roi_idx < sizeByROI.size(); ++roi_idx) {
    unsigned int roi = sizeByROI[roi_idx].second;
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(roi); // the label number is the connected component number - not the one label as mask
    int label = labelObject->GetLabel();

    // do this for the left side
    int targetLocationStart[2];
    targetLocationStart[0] = 0;
    targetLocationStart[1] = roi_idx * base_image_sizeLH;
    int targetLocationEnd[2];
    targetLocationEnd[0] = targetLocationStart[0] + base_image_sizeLW;
    targetLocationEnd[1] = targetLocationStart[1] + base_image_sizeLH;

    RegionType::SizeType targetSize;
    targetSize[0] = targetLocationEnd[0] - targetLocationStart[0];
    targetSize[1] = targetLocationEnd[1] - targetLocationStart[1];
    RegionType::IndexType targetIndex;
    targetIndex[0] = targetLocationStart[0];
    targetIndex[1] = targetLocationStart[1];

    // compute the bounding box for this object in all three orientations
    std::vector<double> boundingBox(6); // three coordinates and one label value (in physical coordinates)
    using PT = typename ImageType3D::PointType;
    PT floatIndexA; // a single point's coordinate (in world coordinates)
    itk::Index<3U> index;
    for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
      index = labelObject->GetIndex(pixelId);

      // get the position in floating point from the index
      image->TransformIndexToPhysicalPoint(index, floatIndexA);
      if (pixelId == 0) {
        boundingBox[0] = floatIndexA[0];
        boundingBox[1] = floatIndexA[1];
        boundingBox[2] = floatIndexA[2];
        boundingBox[3] = floatIndexA[0];
        boundingBox[4] = floatIndexA[1];
        boundingBox[5] = floatIndexA[2];
      }
      if (boundingBox[0] > floatIndexA[0])
        boundingBox[0] = floatIndexA[0];
      if (boundingBox[1] > floatIndexA[1])
        boundingBox[1] = floatIndexA[1];
      if (boundingBox[2] > floatIndexA[2])
        boundingBox[2] = floatIndexA[2];
      if (boundingBox[3] < floatIndexA[0])
        boundingBox[3] = floatIndexA[0];
      if (boundingBox[4] < floatIndexA[1])
        boundingBox[4] = floatIndexA[1];
      if (boundingBox[5] < floatIndexA[2])
        boundingBox[5] = floatIndexA[2];
    }
    // now make the bounding box bigger
    std::vector<double> biggerBB(6);
    // A bounding box can be too small if the label is tiny, we should enlarge the 
    // bounding box suffiently to reach a minimum of maybe 30% of the space in the volume.
    float enlarge[3];
    enlarge[0] = (boundingBox[3]-boundingBox[0]);
    enlarge[1] = (boundingBox[4]-boundingBox[1]);
    enlarge[2] = (boundingBox[5]-boundingBox[2]);
    // if total volume still too small make it even larger, boundingBox size is in mm
    if (enlarge[0]*3 < 100)
      enlarge[0] = 100 - enlarge[0];
    if (enlarge[1]*3 < 100)
      enlarge[1] = 100 - enlarge[1];
    if (enlarge[2]*3 < 100)
      enlarge[2] = 100 - enlarge[2];


    biggerBB[0] = boundingBox[0] - enlarge[0];
    biggerBB[1] = boundingBox[1] - enlarge[1];
    biggerBB[2] = boundingBox[2] - enlarge[2];
    biggerBB[3] = boundingBox[3] + enlarge[0];
    biggerBB[4] = boundingBox[4] + enlarge[1];
    biggerBB[5] = boundingBox[5] + enlarge[2];

    for (int counter = 0; counter < 6; counter++)
      boundingBox[counter] = biggerBB[counter];

    // by definition our label is in the center of the picture
    returns.pos.push_back(std::array<int, 2>{(int)floor(targetLocationStart[0] + (targetLocationEnd[0]-targetLocationStart[0])/2.0), (int)floor(targetLocationStart[1] + (targetLocationEnd[1]-targetLocationStart[1])/2.0)});
    returns.text.push_back(std::to_string(label));

    // do the next steps 3 times, start with left image
    {
      RegionType targetRegion(targetIndex, targetSize); // the region we want to fill in the output key image
      itk::ImageRegionIteratorWithIndex<CImageType> targetRGBIterator(keyImage, targetRegion);

      // finished computing the maximum enclosing bounding box for this region of interest
      // do we need that? We could zoom in, but that is costly... with a zoom factor. 
      // We would want to see the region of interest and double the space around the region.
      targetRGBIterator.GoToBegin(); // 2D Volume of curved slice
      while (!targetRGBIterator.IsAtEnd() ) {
        CPixelType value = targetRGBIterator.Value(); // we ignore this value, will be filled in with correct one at the end
        CImageType::IndexType idx = targetRGBIterator.GetIndex();
        // the current output pixel (idx) in percentage of output image (index coordinates)
        float idx_percent[2];
        idx_percent[0] = ((float)idx[0] - (float)targetLocationStart[0])/((float)targetLocationEnd[0]-(float)targetLocationStart[0]);
        idx_percent[1] = ((float)idx[1] - (float)targetLocationStart[1])/((float)targetLocationEnd[1]-(float)targetLocationStart[1]);
        // the above should be between 0.0 and 1.0, test here
        if (idx_percent[0] < 0.0 || idx_percent[0]>1.0) {
          fprintf(stderr, "Error: idx_percent is not 0..1 but %f. Should not happen\n", idx_percent[0]);
        }
        if (idx_percent[1] < 0.0 || idx_percent[1]>1.0) {
          fprintf(stderr, "Error: idx_percent is not 0..1 but %f. Should not happen\n", idx_percent[1]);
        }

        itk::ContinuousIndex<double, 3> pixel;
        itk::ContinuousIndex<double, 3> floatIndexA;
        // use x as fastest running index and y as second fast running, third dimension is at mid-point
        // TODO: scaling of voxel size is not taken care of. We need to move differently based on pixel size
        //auto voxelSize = image->GetSpacing();
        float start0 = boundingBox[0];
        float start1 = boundingBox[1];
        float start2 = boundingBox[2];
        float sx = (boundingBox[3]-boundingBox[0]);
        float sy = (boundingBox[4]-boundingBox[1]);
        float sz = (boundingBox[5]-boundingBox[2]);
        if (sx > sy) {
          // adjust the other bounding box to the same size, and center it
          float midy = start1 + sy/2.0;
          start1 = midy - (sx/2.0);
          sy = sx;
        } else {
          // adjust the other bounding box to the same size, and center it
          float midx = start0 + sx/2.0;
          start0 = midx - (sy/2.0);
          sx = sy;
        }
        // aspect ratio of output image is 1, so bounding box needs to be adjusted in the smaller dimension to have
        // same aspect ratio based on larger dimension sx or sy
        pixel[0] = (float)start0 + idx_percent[0]*sx; // in pixel coordinates of image, based on boundingBox of one dimension
        pixel[1] = (float)start1 + idx_percent[1]*sy;
        pixel[2] = (float)start2 + sz/2.0f; // middle

        image->TransformPhysicalPointToContinuousIndex(pixel, floatIndexA);

        // std::cout << "Value at 1.3: " << interpolator->EvaluateAtContinuousIndex(pixel) << std::endl;
        float grayValue = 0;
        float redValue = 0;
        float blueValue = 0;
        float greenValue = 0;
        if (interpolator->IsInsideBuffer(floatIndexA)) {
          grayValue = interpolator->EvaluateAtContinuousIndex(floatIndexA);
          redValue = interpolatorRed->EvaluateAtContinuousIndex(floatIndexA);
          greenValue = interpolatorGreen->EvaluateAtContinuousIndex(floatIndexA);
          blueValue = interpolatorBlue->EvaluateAtContinuousIndex(floatIndexA);
        } else {
          //fprintf(stdout, "OUTSIDE region element at location %f %f %f\n", pixel[0], pixel[1], pixel[2]);
          //fflush(stdout);
        }

        float scaledGrayValue = (grayValue - t1) / (t2-t1);

        float red = f * scaledGrayValue + (1 - f) * redValue;
        float green = f * scaledGrayValue + (1 - f) * greenValue;
        float blue = f * scaledGrayValue + (1 - f) * blueValue;
        //fprintf(stdout, "before clipping red, green blue: %f %f %f\n", red, green, blue);
        red = std::min<float>(1, std::max<float>(0,red));
        green = std::min<float>(1, std::max<float>(0,green));
        blue = std::min<float>(1, std::max<float>(0,blue));
        //fprintf(stdout, "red, green blue: %f %f %f\n", red, green, blue);
        //fflush(stdout);
        value.SetRed((int)(red * 255));
        value.SetGreen((int)(green * 255));
        value.SetBlue((int)(blue * 255));
        targetRGBIterator.Set(value);

        ++targetRGBIterator;
      }
    } // Left image done
    { // Right top image (flip 180)

      targetLocationStart[0] = base_image_sizeLW;
      targetLocationStart[1] = roi_idx * base_image_sizeLH;
      targetLocationEnd[0] = targetLocationStart[0] + base_image_sizeRW;
      targetLocationEnd[1] = targetLocationStart[1] + base_image_sizeRH;

      RegionType::SizeType targetSize;
      targetSize[0] = targetLocationEnd[0] - targetLocationStart[0];
      targetSize[1] = targetLocationEnd[1] - targetLocationStart[1];
      RegionType::IndexType targetIndex;
      targetIndex[0] = targetLocationStart[0];
      targetIndex[1] = targetLocationStart[1];

      RegionType targetRegion(targetIndex, targetSize); // the region we want to fill in the output key image
      itk::ImageRegionIteratorWithIndex<CImageType> targetRGBIterator(keyImage, targetRegion);

      // finished computing the maximum inclosing bounding box for this region of interest
      // do we need that? We could zoom in, but that is costly... with a zoom factor. 
      // We would want to see the region of interest and double the space around the region.
      targetRGBIterator.GoToBegin(); // 2D Volume of curved slice
      while (!targetRGBIterator.IsAtEnd() ) {
        CPixelType value = targetRGBIterator.Value(); // we ignore this value, will be filled in with correct one at the end
        CImageType::IndexType idx = targetRGBIterator.GetIndex();
        // the current output pixel (idx) in percentage of output image (index coordinates)
        float idx_percent[2];
        idx_percent[0] = ((float)idx[0] - (float)targetLocationStart[0])/((float)targetLocationEnd[0]-(float)targetLocationStart[0]);
        idx_percent[1] = ((float)idx[1] - (float)targetLocationStart[1])/((float)targetLocationEnd[1]-(float)targetLocationStart[1]);
        idx_percent[1] = 1.0f - idx_percent[1];
        // the above should be between 0.0 and 1.0, test here
        if (idx_percent[0] < 0.0 || idx_percent[0]>1.0) {
          fprintf(stderr, "Error: idx_percent is not 0..1 but %f. Should not happen\n", idx_percent[0]);
        }
        if (idx_percent[1] < 0.0 || idx_percent[1]>1.0) {
          fprintf(stderr, "Error: idx_percent is not 0..1 but %f. Should not happen\n", idx_percent[1]);
        }

        itk::ContinuousIndex<double, 3> pixel;
        itk::ContinuousIndex<double, 3> floatIndexA;

        float start0 = boundingBox[0];
        float start1 = boundingBox[1];
        float start2 = boundingBox[2];
        float sx = (boundingBox[3]-boundingBox[0]);
        float sy = (boundingBox[4]-boundingBox[1]);
        float sz = (boundingBox[5]-boundingBox[2]);
        if (sz > sx) {
          // adjust the other bounding box to the same size, and center it
          float midx = start0 + sx/2.0;
          start0 = midx - (sz/2.0);
          sx = sz;
        } else {
          // adjust the other bounding box to the same size, and center it
          float midz = start2 + sz/2.0;
          start2 = midz - (sx/2.0);
          sz = sx;
        }
        // aspect ratio of output image is 1, so bounding box needs to be adjusted in the smaller dimension to have
        // same aspect ratio based on larger dimension sx or sy
        pixel[0] = (float)start0 + idx_percent[0]*sx; // in pixel coordinates of image, based on boundingBox of one dimension
        pixel[1] = (float)start1 + sy/2.0;
        pixel[2] = (float)start2 + idx_percent[1]*sz; // middle

        // use x as fastest running index and y as second fast running, third dimension is at mid-point
        //pixel[0] = (float)boundingBox[0] + idx_percent[0]*(boundingBox[3]-boundingBox[0]); // in pixel coordinates of image, based on boundingBox of one dimension
        //pixel[1] = (float)boundingBox[1] + (boundingBox[4]-boundingBox[1])/2.0f; // middle
        //pixel[2] = (float)boundingBox[2] + idx_percent[1]*(boundingBox[5]-boundingBox[2]);

        image->TransformPhysicalPointToContinuousIndex(pixel, floatIndexA);

        // std::cout << "Value at 1.3: " << interpolator->EvaluateAtContinuousIndex(pixel) << std::endl;
        float grayValue = 0;
        float redValue = 0;
        float blueValue = 0;
        float greenValue = 0;
        if (interpolator->IsInsideBuffer(floatIndexA)) {
          grayValue = interpolator->EvaluateAtContinuousIndex(floatIndexA);
          redValue = interpolatorRed->EvaluateAtContinuousIndex(floatIndexA);
          greenValue = interpolatorGreen->EvaluateAtContinuousIndex(floatIndexA);
          blueValue = interpolatorBlue->EvaluateAtContinuousIndex(floatIndexA);
        } else {
          //fprintf(stdout, "OUTSIDE region element at location %f %f %f\n", pixel[0], pixel[1], pixel[2]);
          //fflush(stdout);
        }

        float scaledGrayValue = (grayValue - t1) / (t2-t1);

        float red = f * scaledGrayValue + (1 - f) * redValue;
        float green = f * scaledGrayValue + (1 - f) * greenValue;
        float blue = f * scaledGrayValue + (1 - f) * blueValue;
        //fprintf(stdout, "before clipping red, green blue: %f %f %f\n", red, green, blue);
        red = std::min<float>(1, std::max<float>(0,red));
        green = std::min<float>(1, std::max<float>(0,green));
        blue = std::min<float>(1, std::max<float>(0,blue));
        //fprintf(stdout, "red, green blue: %f %f %f\n", red, green, blue);
        //fflush(stdout);
        value.SetRed((int)(red * 255));
        value.SetGreen((int)(green * 255));
        value.SetBlue((int)(blue * 255));
        targetRGBIterator.Set(value);

        ++targetRGBIterator;
      }
    } // Right top image done
    { // Right bottom image

      targetLocationStart[0] = base_image_sizeLW;
      targetLocationStart[1] = roi_idx * base_image_sizeLH + base_image_sizeRH;
      targetLocationEnd[0] = targetLocationStart[0] + base_image_sizeRW;
      targetLocationEnd[1] = targetLocationStart[1] + base_image_sizeRH;

      RegionType::SizeType targetSize;
      targetSize[0] = targetLocationEnd[0] - targetLocationStart[0];
      targetSize[1] = targetLocationEnd[1] - targetLocationStart[1];
      RegionType::IndexType targetIndex;
      targetIndex[0] = targetLocationStart[0];
      targetIndex[1] = targetLocationStart[1];

      RegionType targetRegion(targetIndex, targetSize); // the region we want to fill in the output key image
      itk::ImageRegionIteratorWithIndex<CImageType> targetRGBIterator(keyImage, targetRegion);

      // finished computing the maximum inclosing bounding box for this region of interest
      // do we need that? We could zoom in, but that is costly... with a zoom factor. 
      // We would want to see the region of interest and double the space around the region.
      targetRGBIterator.GoToBegin(); // 2D Volume of curved slice
      while (!targetRGBIterator.IsAtEnd() ) {
        CPixelType value = targetRGBIterator.Value(); // we ignore this value, will be filled in with correct one at the end
        CImageType::IndexType idx = targetRGBIterator.GetIndex();
        // the current output pixel (idx) in percentage of output image (index coordinates)
        float idx_percent[2];
        idx_percent[0] = ((float)idx[0] - (float)targetLocationStart[0])/((float)targetLocationEnd[0]-(float)targetLocationStart[0]);
        idx_percent[1] = ((float)idx[1] - (float)targetLocationStart[1])/((float)targetLocationEnd[1]-(float)targetLocationStart[1]);
        // the above should be between 0.0 and 1.0, test here
        idx_percent[1] = 1.0f - idx_percent[1];
        if (idx_percent[0] < 0.0 || idx_percent[0]>1.0) {
          fprintf(stderr, "Error: idx_percent is not 0..1 but %f. Should not happen\n", idx_percent[0]);
        }
        if (idx_percent[1] < 0.0 || idx_percent[1]>1.0) {
          fprintf(stderr, "Error: idx_percent is not 0..1 but %f. Should not happen\n", idx_percent[1]);
        }

        itk::ContinuousIndex<double, 3> pixel;
        itk::ContinuousIndex<double, 3> floatIndexA;

        float start0 = boundingBox[0];
        float start1 = boundingBox[1];
        float start2 = boundingBox[2];
        float sx = (boundingBox[3]-boundingBox[0]);
        float sy = (boundingBox[4]-boundingBox[1]);
        float sz = (boundingBox[5]-boundingBox[2]);
        if (sz > sy) {
          // adjust the other bounding box to the same size, and center it
          float midy = start1 + sy/2.0;
          start1 = midy - (sz/2.0);
          sy = sz;
        } else {
          // adjust the other bounding box to the same size, and center it
          float midz = start2 + sz/2.0;
          start2 = midz - (sy/2.0);
          sz = sy;
        }
        // aspect ratio of output image is 1, so bounding box needs to be adjusted in the smaller dimension to have
        // same aspect ratio based on larger dimension sx or sy
        pixel[0] = (float)start0 + sx/2.0;; // in pixel coordinates of image, based on boundingBox of one dimension
        pixel[1] = (float)start1 + idx_percent[0]*sy;
        pixel[2] = (float)start2 + idx_percent[1]*sz; // middle

        // use x as fastest running index and y as second fast running, third dimension is at mid-point
        //pixel[0] = (float)boundingBox[0] + (boundingBox[3]-boundingBox[0])/2.0f; // middle
        //pixel[1] = (float)boundingBox[1] + idx_percent[0]*(boundingBox[4]-boundingBox[1]); // in pixel coordinates of image, based on boundingBox of one dimension
        //pixel[2] = (float)boundingBox[2] + idx_percent[1]*(boundingBox[5]-boundingBox[2]);

        image->TransformPhysicalPointToContinuousIndex(pixel, floatIndexA);

        // std::cout << "Value at 1.3: " << interpolator->EvaluateAtContinuousIndex(pixel) << std::endl;
        float grayValue = 0;
        float redValue = 0;
        float blueValue = 0;
        float greenValue = 0;
        if (interpolator->IsInsideBuffer(floatIndexA)) {
          grayValue = interpolator->EvaluateAtContinuousIndex(floatIndexA);
          redValue = interpolatorRed->EvaluateAtContinuousIndex(floatIndexA);
          greenValue = interpolatorGreen->EvaluateAtContinuousIndex(floatIndexA);
          blueValue = interpolatorBlue->EvaluateAtContinuousIndex(floatIndexA);
        } else {
          //fprintf(stdout, "OUTSIDE region element at location %f %f %f\n", pixel[0], pixel[1], pixel[2]);
          //fflush(stdout);
        }

        float scaledGrayValue = (grayValue - t1) / (t2-t1);

        float red = f * scaledGrayValue + (1 - f) * redValue;
        float green = f * scaledGrayValue + (1 - f) * greenValue;
        float blue = f * scaledGrayValue + (1 - f) * blueValue;
        //fprintf(stdout, "before clipping red, green blue: %f %f %f\n", red, green, blue);
        red = std::min<float>(1, std::max<float>(0,red));
        green = std::min<float>(1, std::max<float>(0,green));
        blue = std::min<float>(1, std::max<float>(0,blue));
        //fprintf(stdout, "red, green blue: %f %f %f\n", red, green, blue);
        //fflush(stdout);
        value.SetRed((int)(red * 255));
        value.SetGreen((int)(green * 255));
        value.SetBlue((int)(blue * 255));
        targetRGBIterator.Set(value);

        ++targetRGBIterator;
      }
    } // Right bottom image done

  } 

  return returns;
}


// generate a key image based on the input image and a mask (fused with names)
// test: 
//     ./imageAndMask2Report data/ror_trigger_run_Wednesday_980595789/ror_trigger_run_Wednesday_980595789/input data/ror_trigger_run_Wednesday_980595789/ror_trigger_run_Wednesday_980595789_output/labels/508bc8c54546f0c3383f4325ec6fa70e310328932af7bffcf812079391445.1/ /tmp/bla -u | less
generateImageReturn generateKeyImage(ImageType3D::Pointer image, LabelMapType *labelMap, std::vector<int> resolution, float lowerT, float upperT) {
  if (verbose) {
    fprintf(stdout, "Start generating a key image...\n");
  }
  std::vector<std::vector<float>> labelColors2 = {{0, 0, 0}, {166,206,227}, {31,120,180}, {178,223,138}, {51,160,44}, {251,154,153}, {227,26,28}, {253,191,111}, {255,127,0}, {202,178,214}, {106,61,154}, {255,255,153}, {177,89,40}};

  generateImageReturn returns; // store keyImage and the location and text that should be presented ontop of the image

  // create a new RGB image
  CImageType::Pointer keyImage = CImageType::New();
  returns.keyImage = keyImage;

  using RegionType = itk::ImageRegion<2>;
  RegionType::SizeType size;
  size[0] = resolution[0]; // 512
  size[1] = resolution[1]; // 512

  RegionType::IndexType index;
  index.Fill(0);

  RegionType region(index, size);
  keyImage->SetRegions(region);

  keyImage->Allocate();
  keyImage->FillBuffer(itk::NumericTraits<CPixelType>::Zero);

  //
  // for each label compute the center of mass for a spline function
  //
  std::vector< std::array<double, 4> > centers; // three coordinates and one label value (in physical coordinates)
  std::array<double, 3> minPos, maxPos;
  // keep a map of all pixel coordinates with label
  //std::map< std::string, int > coord2Label;

  for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n) {
    // these are all in order
    returns.roi_order.push_back(n);
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n); // the label number is the connected component number - not the one label as mask
    int label = labelObject->GetLabel();

    using PT = typename ImageType3D::PointType;
    PT floatIndexA; // a single points coordinates in world coordinates
    PT centerHere; // we accumulate positions here to compute the center of mass for this object in 3D (world coordinates)
    centerHere[0] = 0.0f;
    centerHere[1] = 0.0f;
    centerHere[2] = 0.0f;
    itk::Index<3U> index;
    for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
      index = labelObject->GetIndex(pixelId);
      // if we compute the key like this we can get the label for it,
      // not being in the map means that we have a background pixel
      //std::string key = std::to_string(index[0]) + "_" + std::to_string(index[1]) + "_" + std::to_string(index[2]);
      //coord2Label.insert( std::make_pair(key, label) );

      // get the position in floating point from the index
      image->TransformIndexToPhysicalPoint(index, floatIndexA);
      centerHere[0] += floatIndexA[0];
      centerHere[1] += floatIndexA[1];
      centerHere[2] += floatIndexA[2];
    }
    centerHere[0] /= labelObject->Size(); // center of mass
    centerHere[1] /= labelObject->Size();
    centerHere[2] /= labelObject->Size();
    centers.push_back(std::array<double, 4>{centerHere[0], centerHere[1], centerHere[2], (double)label}); // keep track of the label (connected component generated int value)
    if (n == 0) { // init
      minPos[0] = centerHere[0]; minPos[1] = centerHere[1]; minPos[2] = centerHere[2];
      maxPos[0] = centerHere[0]; maxPos[1] = centerHere[1]; maxPos[2] = centerHere[2];
    } else {
      for (int i = 0; i < 3; i++) {
        if (centerHere[i] < minPos[i]) {
          minPos[i] = centerHere[i];
        }
        if (centerHere[i] > maxPos[i]) {
          maxPos[i] = centerHere[i];
        }
      } 
    }
  }
  if (centers.size() == 0) {
    // nothing to do, return here
    if (verbose)
      fprintf(stdout, "nothing to do here, no centers found\n");
    return returns;
  }
  // compute the extend in all three dimensions and the index of the longest axis
  std::vector<double> dists{(maxPos[0]-minPos[0]),(maxPos[1]-minPos[1]),(maxPos[2]-minPos[2])};
  int directionLongestAxis = std::max_element(dists.begin(), dists.end()) - dists.begin();
  
  if (verbose)
    fprintf(stdout, "direction longest axis is: %d\n", directionLongestAxis);

  /*if (verbose)
    for (int i = 0; i < centers.size(); i++) {
      fprintf(stdout, "centers before sort %d is : %f %f %f\n", i, centers[i][0], centers[i][1], centers[i][2]);
    }*/

  // we want to interpolate between the points in the right order along the directionLongestAxis
  // TODO: find out what the correct axis is based on some DICOM tags
  std::sort(centers.begin(), centers.end(), [directionLongestAxis](const std::array<double, 4> a, const std::array<double, 4> b) {
    return a.at(directionLongestAxis) > b.at(directionLongestAxis);
  });

  /*if (verbose)
    for (int i = 0; i < centers.size(); i++) {
      fprintf(stdout, "centers after sort %d is : %f %f %f\n", i, centers[i][0], centers[i][1], centers[i][2]);
    }*/

  // add two points at the beginning and at the end to continue to the center curve
  auto one = centers[0];
  // we can guarantee that we have one region here
  if (centers.size() > 1) {
    auto two = centers[1];
    std::array<double, 4> newBeginning = std::array<double, 4>{ one[0] + (one[0] - two[0]), one[1] + (one[1] - two[1]), one[2] + (one[2] - two[2]), -1 };
    centers.insert(centers.begin(), newBeginning);

    one = centers[centers.size()-1];
    two = centers[centers.size()-2];
    std::array<double, 4> newEnding = std::array<double, 4>{ one[0] + (one[0] - two[0]), one[1] + (one[1] - two[1]), one[2] + (one[2] - two[2]), (double)(centers.size()+1) };
    centers.push_back(newEnding);
  } else {
    auto two = std::array<double, 4>{ centers[0][0], centers[0][1], centers[0][2], -1 };
    two[directionLongestAxis] += 10.0;
    std::array<double, 4> newBeginning = std::array<double, 4>{ one[0] + (one[0] - two[0]), one[1] + (one[1] - two[1]), one[2] + (one[2] - two[2]), -1 };
    centers.insert(centers.begin(), newBeginning);

    one = centers[centers.size()-1];
    two = centers[centers.size()-2];
    std::array<double, 4> newEnding = std::array<double, 4>{ one[0] + (one[0] - two[0]), one[1] + (one[1] - two[1]), one[2] + (one[2] - two[2]), (double)(centers.size()+1) };
    centers.push_back(newEnding);
  }

  // create an open spline
  // Warning: we need at least 4 objects here
  // If we have less than 4 but more than 2 regions we can use a linear function from start to end and hope for the best


  // we can evaluate the spline
  // TODO: use a linear extension in the direction of cr.prime() (tangent vector)
  std::vector< std::vector< double > > interpolatedCenterLocations; // should be resolution[1] many, in physical space
  if (centers.size() > 4) {
    std::vector< std::array<double, 3> > positions; // make a copy of the first 3 dimensions
    for (int p = 0; p < centers.size(); p++) {
      positions.push_back( std::array<double, 3>{centers[p][0], centers[p][1], centers[p][2]});
    }
    boost::math::catmull_rom<std::array<double, 3>> cr(std::move(positions));
    double dt = cr.max_parameter()/((float)resolution[1]-1.0);
    for (int p = 0; p < resolution[1]; p++) {
      float s = (double)p * dt;
      auto point = cr( (double)p * dt );
      interpolatedCenterLocations.push_back(std::vector<double>{point[0], point[1], point[2]});
      //fprintf(stdout, "found a point at %d %f here: %f %f %f\n", p, s, point[0], point[1], point[2]);
    }
    //fflush(stdout);
  } else if (centers.size() > 1) {
    for (int p = 0; p < resolution[1]; p++) {
      std::vector<double> point;
      point.push_back(centers[0][0] + p*(centers[centers.size()-1][0] - centers[0][0])/((double)resolution[1]-1.0)  );
      point.push_back(centers[0][1] + p*(centers[centers.size()-1][1] - centers[0][1])/((double)resolution[1]-1.0)  );
      point.push_back(centers[0][2] + p*(centers[centers.size()-1][2] - centers[0][2])/((double)resolution[1]-1.0)  );
      interpolatedCenterLocations.push_back(std::vector<double>{point[0], point[1], point[2]});
    }
  } else {
    fprintf(stderr, "Error: no key image for less than 2 points!\n");
    // TODO: do a 3 axis view
  }

  // now sample the image, instead of normal and binormal use one of the other two dimensions
  // sample dimension is directionLongestAxis+1 % 3
  std::vector<double> sampleDirection{0,0,0};
  int sampleDimension = (directionLongestAxis-1) % 3;
  if (sampleDimension < 0 || sampleDimension > 2)
    sampleDimension = 0; // fallback
  sampleDirection[sampleDimension] = 1.0;
  double stepSize = sqrtf( /*(interpolatedCenterLocations[1][0] - interpolatedCenterLocations[2][0]) * (interpolatedCenterLocations[1][0] - interpolatedCenterLocations[2][0]) +  
                           (interpolatedCenterLocations[1][1] - interpolatedCenterLocations[2][1]) * (interpolatedCenterLocations[1][1] - interpolatedCenterLocations[2][1]) + */ 
                           (interpolatedCenterLocations[1][directionLongestAxis] - interpolatedCenterLocations[2][directionLongestAxis]) * (interpolatedCenterLocations[1][directionLongestAxis] - interpolatedCenterLocations[2][directionLongestAxis]) );
  if (verbose) {
    fprintf(stdout, "sample direction: %f %f %f, stepsize: %f\n", sampleDirection[0], sampleDirection[1], sampleDirection[2], stepSize);
    fflush(stdout);
  }
  // now sample the output image using the coordinates system we established above
  using IteratorTypeImage = itk::ImageRegionIteratorWithIndex< ImageType3D >;
  // we can use coord2Label to lookup the mask values, maybe color is missing?
  // we need to go through each pixel of the output image

  RegionType outputRegion = keyImage->GetLargestPossibleRegion();
  itk::ImageRegionIteratorWithIndex<CImageType> outputRGBIterator(keyImage, outputRegion);

  //
  // compute optimal window level for whole image, we will use the values for the curved slice
  //
  using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<ImageType3D>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(image);
  imageCalculatorFilter->Compute();
  int minGray = imageCalculatorFilter->GetMinimum();
  int maxGray = imageCalculatorFilter->GetMaximum();

  using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImageType3D>;
  HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(image);
  int histogramSize = 1024;
  histogramGenerator->SetNumberOfBins(histogramSize);
  histogramGenerator->SetHistogramMin(minGray);
  histogramGenerator->SetHistogramMax(maxGray);
  histogramGenerator->Compute();
  using HistogramType = HistogramGeneratorType::HistogramType;
  const HistogramType *histogram = histogramGenerator->GetOutput();
  // set in calling function
  //double lowerT = 0.01;
  //double upperT = 0.999;
  double t1 = -1;
  double t2 = -1;
  double sum = 0;
  double total = 0;
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    total += histogram->GetFrequency(bin, 0);
  }
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    double f = histogram->GetFrequency(bin, 0) / total;
    // fprintf(stdout, "bin %d, value is %f\n", bin, f);
    sum += f;
    if (t1 == -1 && sum > lowerT) {
      t1 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
    }
    if (t2 == -1 && sum > upperT) {
      t2 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
      break;
    }
  }
  //if (verbose) {
  //  fprintf(stdout, "calculated best threshold low: %f, high: %f\n", t1, t2);
  //}

  // we need three 3D images for the red green and blue channel
  // so we can blurr them before using them in the output image
  typedef float FPixelType;
  // typedef itk::Image<FPixelType, 2> FloatImageType;
  using FloatImageType = itk::Image<FPixelType, 3>;
  FloatImageType::Pointer red_channel = FloatImageType::New();
  FloatImageType::Pointer green_channel = FloatImageType::New();
  FloatImageType::Pointer blue_channel = FloatImageType::New();
  ImageType3D::RegionType fusedRegion = image->GetLargestPossibleRegion();
  red_channel->SetRegions(fusedRegion);
  red_channel->Allocate();
  red_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  red_channel->SetOrigin(image->GetOrigin());
  red_channel->SetSpacing(image->GetSpacing());
  red_channel->SetDirection(image->GetDirection());
  green_channel->SetRegions(fusedRegion);
  green_channel->Allocate();
  green_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  green_channel->SetOrigin(image->GetOrigin());
  green_channel->SetSpacing(image->GetSpacing());
  green_channel->SetDirection(image->GetDirection());
  blue_channel->SetRegions(fusedRegion);
  blue_channel->Allocate();
  blue_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  blue_channel->SetOrigin(image->GetOrigin());
  blue_channel->SetSpacing(image->GetSpacing());
  blue_channel->SetDirection(image->GetDirection());
//  itk::ImageRegionIterator<FloatImageType> redIterator(red_channel, fusedRegion);
//  itk::ImageRegionIterator<FloatImageType> greenIterator(green_channel, fusedRegion);
//  itk::ImageRegionIterator<FloatImageType> blueIterator(blue_channel, fusedRegion);

  for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); n++) {
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n); // the label number is the connected component number - not the one label as mask
    int label = labelObject->GetLabel(); // it might be that all regions have the same label and only n is a good choice for the color

    // color is
    std::vector<float> col = labelColors2[ (label % (labelColors2.size()-1)) +1 ];
    itk::Index<3U> index;
    for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
      index = labelObject->GetIndex(pixelId);
      // set this position in all three color images, values are between 0 and 1
      red_channel->SetPixel(index, col[0]/255.0); 
      green_channel->SetPixel(index, col[1]/255.0);
      blue_channel->SetPixel(index, col[2]/255.0);
    }
  }

  // smooth all three channels independently from each other
  using GFilterType = itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>;
  auto gaussFilterR = GFilterType::New();
  gaussFilterR->SetInput(red_channel);
  gaussFilterR->SetVariance(1.5f);
  gaussFilterR->Update();
  FloatImageType::Pointer smoothRed = gaussFilterR->GetOutput();

  auto gaussFilterG = GFilterType::New();
  gaussFilterG->SetInput(green_channel);
  gaussFilterG->SetVariance(1.5f);
  gaussFilterG->Update();
  FloatImageType::Pointer smoothGreen = gaussFilterG->GetOutput();

  auto gaussFilterB = GFilterType::New();
  gaussFilterB->SetInput(blue_channel);
  gaussFilterB->SetVariance(1.5f);
  gaussFilterB->Update();
  FloatImageType::Pointer smoothBlue = gaussFilterB->GetOutput();

  itk::ImageRegionIterator<FloatImageType> redSIterator(smoothRed, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> greenSIterator(smoothGreen, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> blueSIterator(smoothBlue, fusedRegion);

//  redSIterator.GoToBegin(); // 3D Volume iterators
//  greenSIterator.GoToBegin();
//  blueSIterator.GoToBegin();

  // we want to do a tri-linear lookup in the input image
  itk::LinearInterpolateImageFunction<ImageType3D, double>::Pointer interpolator =
    itk::LinearInterpolateImageFunction<ImageType3D, double>::New();
  interpolator->SetInputImage(image);

  itk::LinearInterpolateImageFunction<FloatImageType, double>::Pointer interpolatorRed =
    itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolatorRed->SetInputImage(smoothRed);

  itk::LinearInterpolateImageFunction<FloatImageType, double>::Pointer interpolatorGreen =
    itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolatorGreen->SetInputImage(smoothGreen);

  itk::LinearInterpolateImageFunction<FloatImageType, double>::Pointer interpolatorBlue =
    itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolatorBlue->SetInputImage(smoothBlue);

  float f = 0.7;
  std::vector<int> centerWinner;
  // we ignore the first and last element (added but not really centers of segmented spines)
  for (int p = 1; p < centers.size()-1; p++) {
    // for each of the centers, find the index of the closest point in interpolatedCenterLocations
    float dist = 100000;
    int winner = 0;
    for (int q = 0; q < interpolatedCenterLocations.size(); q++) {
      float current = sqrtf( (centers[p][0] - interpolatedCenterLocations[q][0])*(centers[p][0] - interpolatedCenterLocations[q][0]) + 
                             (centers[p][1] - interpolatedCenterLocations[q][1])*(centers[p][1] - interpolatedCenterLocations[q][1]) + 
                             (centers[p][2] - interpolatedCenterLocations[q][2])*(centers[p][2] - interpolatedCenterLocations[q][2]) );
      if (dist > current) {
        //fprintf(stdout, "found a better winner %d with smaller distance: %f\n", q, current);
        dist = current;
        winner = q;
      }
    }
    centerWinner.push_back(winner);
  }

  // TODO: Instead of a straight representation of the curvilinear slice we 
  //       rather like a straightening in one image plane only (L/R).
  //       We need to compute the center axis that is in the middle of all the points, just the mean of all sample points.

  std::vector<double> centralAxis{0,0,0}; // the mid point of all centers, any axis should go throught that point (in some dimension)
  for (int p = 0; p < centers.size(); p++) {
      centralAxis[0] += centers[p][0];
      centralAxis[1] += centers[p][1];
      centralAxis[2] += centers[p][2];
  }
  centralAxis[0] /= centers.size();
  centralAxis[1] /= centers.size();
  centralAxis[2] /= centers.size();

  for (int vert_pos = 0; vert_pos < resolution[1]; vert_pos++) {
    // We need to find the pixel location in outputRGB that matches with the centers to be able
    // to place the texts.
    // The index of the matching point is now in centerWinner
    for (int p = 0; p < centerWinner.size(); p++) {
      if (vert_pos == centerWinner[p]) {
        // add a text to the current location in the image (idx[0], idx[1]) 
        // p is the centers points for which this is the closest id
        // so centers[p][3] is the label for this point
        //fprintf(stdout, "found a center winner at position: %d, %d\n", vert_pos, p);
        //fflush(stdout);
        returns.pos.push_back(std::array<int, 2>{(int)floor(resolution[0]/2.0), vert_pos});
        returns.text.push_back(std::to_string((int)centers[p+1][3]));
        break;
      }
    }
  }

  float scaling = 1.0f; // 3*(image->GetSpacing()[directionLongestAxis ]/(1.0f * image->GetSpacing()[sampleDimension]));
  scaling = 512.0/2.0 * image->GetSpacing()[sampleDimension]/fusedRegion.GetSize()[sampleDimension]; //*image->GetSpacing()[directionLongestAxis];

  outputRGBIterator.GoToBegin(); // 2D Volume of curved slice
  while (!outputRGBIterator.IsAtEnd() ) {
    CPixelType value = outputRGBIterator.Value(); // we ignore this value, will be filled in with correct one at the end
    CImageType::IndexType idx = outputRGBIterator.GetIndex(); // a 2D index
    // compute for this pixel in the image the interpolated values from input and from the label
    std::vector<double> rowCenter = interpolatedCenterLocations[idx[1]];
    // step is idx[0], but 0 is minus half the steps
    int step = idx[0] - floor(resolution[0]/2.0); // -256..256

    // TODO: We would like to keep the shape of the spine and only center along the L/R orientation.
    //       For that we need to change the rowCenter and shift it.

    // the location in input we want to sample for this pixel in the output 2D image, this is in physical space
    // as sampleDirection is zero in two places and 1.0 in another we can multiply to select a shift in one axis only
    std::vector<double> sampleLocation{ rowCenter[0] + (step*scaling) * sampleDirection[0],
                                        rowCenter[1] + (step*scaling) * sampleDirection[1],
                                        rowCenter[2] + (step*scaling) * sampleDirection[2] };


    // sample at this point
    itk::ContinuousIndex<double, 3> pixel;
    itk::ContinuousIndex<double, 3> floatIndexA;
    pixel[0] = sampleLocation[0];
    pixel[1] = sampleLocation[1];
    pixel[2] = sampleLocation[2];

    image->TransformPhysicalPointToContinuousIndex(pixel, floatIndexA);


    // std::cout << "Value at 1.3: " << interpolator->EvaluateAtContinuousIndex(pixel) << std::endl;
    float grayValue = 0;
    float redValue = 0;
    float blueValue = 0;
    float greenValue = 0;
    if (interpolator->IsInsideBuffer(floatIndexA)) {
      grayValue = interpolator->EvaluateAtContinuousIndex(floatIndexA);
      redValue = interpolatorRed->EvaluateAtContinuousIndex(floatIndexA);
      greenValue = interpolatorGreen->EvaluateAtContinuousIndex(floatIndexA);
      blueValue = interpolatorBlue->EvaluateAtContinuousIndex(floatIndexA);
    } else {
      //fprintf(stdout, "OUTSIDE region element at location %f %f %f\n", pixel[0], pixel[1], pixel[2]);
      //fflush(stdout);
    }

    float scaledGrayValue = (grayValue - t1) / (t2-t1);

    float red = f * scaledGrayValue + (1 - f) * redValue;
    float green = f * scaledGrayValue + (1 - f) * greenValue;
    float blue = f * scaledGrayValue + (1 - f) * blueValue;
    //fprintf(stdout, "before clipping red, green blue: %f %f %f\n", red, green, blue);
    red = std::min<float>(1, std::max<float>(0,red));
    green = std::min<float>(1, std::max<float>(0,green));
    blue = std::min<float>(1, std::max<float>(0,blue));
    //fprintf(stdout, "red, green blue: %f %f %f\n", red, green, blue);
    //fflush(stdout);
    value.SetRed((int)(red * 255));
    value.SetGreen((int)(green * 255));
    value.SetBlue((int)(blue * 255));
    outputRGBIterator.Set(value);

    ++outputRGBIterator;
  }

  return returns;
}


void writeSecondaryCapture(MaskImageType2D::Pointer maskFromPolys, std::string filename, std::string p_out, bool uidFixedFlag,
                           std::string newFusedSeriesInstanceUID, std::string newFusedSOPInstanceUID, bool verbose, float lowerT, float upperT) {

  typedef itk::ImageFileReader<ImageType2D> ReaderType;
  ReaderType::Pointer r = ReaderType::New();
  r->SetFileName(filename);

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
  r->SetImageIO(gdcmImageIO);

  try {
    r->Update();
  } catch (itk::ExceptionObject &e) {
    std::cerr << "exception in file reader " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return;
  }

  itk::MetaDataDictionary &dictionarySlice = r->GetOutput()->GetMetaDataDictionary();

  std::string SeriesNumber("");
  std::string studyID(""); // this is the StudyInstanceUID
  std::string AcquisitionNumber("");
  std::string InstanceNumber("");
  std::string frameOfReferenceUID("");
  std::string InstitutionName("");
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|000d", studyID);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0080", InstitutionName);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0011", SeriesNumber);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0012", AcquisitionNumber);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0013", InstanceNumber); // keep that number
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0052", frameOfReferenceUID);

  int newSeriesNumber = atoi(SeriesNumber.c_str()) + 1003;

  ImageType2D::Pointer im2change = r->GetOutput();
  ImageType2D::RegionType region;
  region = im2change->GetBufferedRegion();
  ImageType2D::SizeType size = region.GetSize();

  ImageType2D::PixelContainer *container;
  container = im2change->GetPixelContainer();
  container->SetContainerManageMemory(false);
  // unsigned int bla = sizeof(InputImageType::PixelType);
  ImageType2D::PixelType *buffer2 = container->GetBufferPointer();

  gdcm::PixelFormat pf = gdcm::PixelFormat::UINT8;
  pf.SetSamplesPerPixel(3);
  gdcm::SmartPointer<gdcm::Image> simage = new gdcm::Image;
  gdcm::Image &image = *simage;
  image.SetNumberOfDimensions(2);
  // typedef itk::Image<PixelType, 2> ImageType2D;
  ImageType2D::RegionType inRegion = im2change->GetLargestPossibleRegion();
  image.SetDimension(0, static_cast<unsigned int>(inRegion.GetSize()[0]));
  image.SetDimension(1, static_cast<unsigned int>(inRegion.GetSize()[1]));
  // image.SetDimension(2, m_Dimensions[2] );
  image.SetSpacing(0, im2change->GetSpacing()[0]);
  image.SetSpacing(1, im2change->GetSpacing()[1]);

  image.SetPixelFormat(pf);
  gdcm::PhotometricInterpretation pi = gdcm::PhotometricInterpretation::RGB;
  image.SetPhotometricInterpretation(pi);
  image.SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
  // copy the DICOM tags over from inputImage to image
  gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));

  // create a fused volume using im2change and maskFromPolys
  using CPixelType = itk::RGBPixel<unsigned char>;
  using CImageType = itk::Image<CPixelType, 2>;
  using CWriterType = itk::ImageFileWriter<CImageType>;
  CImageType::Pointer fused = CImageType::New();
  CImageType::RegionType fusedRegion = im2change->GetLargestPossibleRegion();
  fused->SetRegions(fusedRegion);
  fused->Allocate();
  fused->FillBuffer(itk::NumericTraits<CPixelType>::Zero);
  fused->SetOrigin(im2change->GetOrigin());
  fused->SetSpacing(im2change->GetSpacing());
  fused->SetDirection(im2change->GetDirection());

  itk::ImageRegionIterator<ImageType2D> fusedLabelIterator(maskFromPolys, fusedRegion);
  itk::ImageRegionIterator<ImageType2D> inputIterator(im2change, fusedRegion);
  itk::ImageRegionIterator<CImageType> fusedIterator(fused, fusedRegion);

  using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<ImageType2D>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(im2change);
  imageCalculatorFilter->Compute();
  int minGray = imageCalculatorFilter->GetMinimum();
  int maxGray = imageCalculatorFilter->GetMaximum();
  // what are good contrast and brightness values?
  // compute a histogram first
  using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImageType2D>;
  HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(im2change);
  int histogramSize = 1024;
  histogramGenerator->SetNumberOfBins(histogramSize);
  histogramGenerator->SetHistogramMin(minGray);
  histogramGenerator->SetHistogramMax(maxGray);
  histogramGenerator->Compute();
  using HistogramType = HistogramGeneratorType::HistogramType;
  const HistogramType *histogram = histogramGenerator->GetOutput();
  //double lowerT = 0.01;
  //double upperT = 0.999;
  double t1 = -1;
  double t2 = -1;
  double sum = 0;
  double total = 0;
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    total += histogram->GetFrequency(bin, 0);
  }
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    double f = histogram->GetFrequency(bin, 0) / total;
    // fprintf(stdout, "bin %d, value is %f\n", bin, f);
    sum += f;
    if (t1 == -1 && sum > lowerT) {
      t1 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
    }
    if (t2 == -1 && sum > upperT) {
      t2 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
      break;
    }
  }
  if (verbose) {
    fprintf(stdout, "calculated best threshold low: %f, high: %f\n", t1, t2);
  }

  std::vector<std::vector<float>> labelColors = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  std::vector<std::vector<float>> labelColors2 = {{0, 0, 0}, {166,206,227}, {31,120,180}, {178,223,138}, {51,160,44}, {251,154,153}, {227,26,28}, {253,191,111}, {255,127,0}, {202,178,214}, {106,61,154}, {255,255,153}, {177,89,40}};

  // do this in two steps, first compute three label color channels, smooth them and alpha-blend last
  typedef float FPixelType;
  // typedef itk::Image<FPixelType, 2> FloatImageType;
  using FloatImageType = itk::Image<FPixelType, 2>;
  FloatImageType::Pointer red_channel = FloatImageType::New();
  FloatImageType::Pointer green_channel = FloatImageType::New();
  FloatImageType::Pointer blue_channel = FloatImageType::New();
  fusedRegion = im2change->GetLargestPossibleRegion();
  red_channel->SetRegions(fusedRegion);
  red_channel->Allocate();
  red_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  red_channel->SetOrigin(im2change->GetOrigin());
  red_channel->SetSpacing(im2change->GetSpacing());
  red_channel->SetDirection(im2change->GetDirection());
  green_channel->SetRegions(fusedRegion);
  green_channel->Allocate();
  green_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  green_channel->SetOrigin(im2change->GetOrigin());
  green_channel->SetSpacing(im2change->GetSpacing());
  green_channel->SetDirection(im2change->GetDirection());
  blue_channel->SetRegions(fusedRegion);
  blue_channel->Allocate();
  blue_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  blue_channel->SetOrigin(im2change->GetOrigin());
  blue_channel->SetSpacing(im2change->GetSpacing());
  blue_channel->SetDirection(im2change->GetDirection());
  itk::ImageRegionIterator<FloatImageType> redIterator(red_channel, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> greenIterator(green_channel, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> blueIterator(blue_channel, fusedRegion);

  fusedLabelIterator.GoToBegin();
  redIterator.GoToBegin();
  greenIterator.GoToBegin();
  blueIterator.GoToBegin();
  // create a map for label by color 
  std::map<int, int> labelToColor;
  labelToColor.insert(std::make_pair<int, int>(0, 0)); // background is the first label
  while (!fusedLabelIterator.IsAtEnd() && !redIterator.IsAtEnd() && !greenIterator.IsAtEnd() && !blueIterator.IsAtEnd()) {
    // this will crash with many labels (more than in our const array)
    int vvvv = fusedLabelIterator.Get();
    if (labelToColor.find(vvvv) == labelToColor.end()) {
      // get the next color from labelColor2
      if (verbose) {
        fprintf(stdout, " Found a new label %d, set to color: %zd\n", vvvv, (labelToColor.size() % (labelColors2.size()-1))+1 );
        fflush(stdout);
      }
      labelToColor.insert(std::pair<int, int>(vvvv, (labelToColor.size() % (labelColors2.size()-1))+1));
    }
    //if (vvvv < 0 || vvvv > labelColors.size() - 1) {
    //  fprintf(stderr, "Warning: mask label is too large for our colors %d\n", vvvv);
    //}
    std::vector<float> col = labelColors2[labelToColor[vvvv]];
    redIterator.Set(col[0]/255.0); // values are 0..1
    greenIterator.Set(col[1]/255.0);
    blueIterator.Set(col[2]/255.0);
    ++redIterator;
    ++greenIterator;
    ++blueIterator;
    ++fusedLabelIterator;
  }

  // now smooth with Gaussian (each channel independently)
  using GFilterType = itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>;
  auto gaussFilterR = GFilterType::New();
  gaussFilterR->SetInput(red_channel);
  // gaussFilterR->SetMaximumKernelWidth(3);
  gaussFilterR->SetVariance(1.5f);
  gaussFilterR->Update();
  FloatImageType::Pointer smoothRed = gaussFilterR->GetOutput();

  auto gaussFilterG = GFilterType::New();
  gaussFilterG->SetInput(green_channel);
  // gaussFilterG->SetMaximumKernelWidth(3);
  gaussFilterG->SetVariance(1.5f);
  gaussFilterG->Update();
  FloatImageType::Pointer smoothGreen = gaussFilterG->GetOutput();

  auto gaussFilterB = GFilterType::New();
  gaussFilterB->SetInput(blue_channel);
  // gaussFilterB->SetMaximumKernelWidth(3);
  gaussFilterB->SetVariance(1.5f);
  gaussFilterB->Update();
  FloatImageType::Pointer smoothBlue = gaussFilterB->GetOutput();

  itk::ImageRegionIterator<FloatImageType> redSIterator(smoothRed, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> greenSIterator(smoothGreen, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> blueSIterator(smoothBlue, fusedRegion);

  // now use the smaoothed color channels (clamp them between 0 and 1)
  inputIterator.GoToBegin();
  fusedIterator.GoToBegin();
  redSIterator.GoToBegin();
  greenSIterator.GoToBegin();
  blueSIterator.GoToBegin();
  float f = 0.6; // weight of the underlay, at 0.1 only mask is visible
  float red, green, blue;
  while (!inputIterator.IsAtEnd() && !fusedIterator.IsAtEnd() && !redSIterator.IsAtEnd() && !greenSIterator.IsAtEnd() && !blueSIterator.IsAtEnd()) {
    float scaledP = ((float) inputIterator.Get() - t1) / (t2 - t1);
    CPixelType value = fusedIterator.Value();

    red = redSIterator.Get();
    green = greenSIterator.Get();
    blue = blueSIterator.Get();

    // alpha blend
    red = f * scaledP + red * (1 - f);
    green = f * scaledP + green * (1 - f);
    blue = f * scaledP + blue * (1 - f);

    red = std::min<float>(1, std::max<float>(0, red));
    green = std::min<float>(1, std::max<float>(0, green));
    blue = std::min<float>(1, std::max<float>(0, blue));

    value.SetRed((int)(red * 255));
    value.SetGreen((int)(green * 255));
    value.SetBlue((int)(blue * 255));
    fusedIterator.Set(value);

    ++inputIterator; 
    ++fusedIterator;
    ++redSIterator;
    ++greenSIterator;
    ++blueSIterator;    
  }

  // now change the DICOM tags for the series and save it again
  // itk::MetaDataDictionary &dictionarySlice = r->GetOutput()->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000d", studyID); // StudyInstanceUID
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0010", studyID); // StudyID  
  //  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000E", newFusedSeriesInstanceUID); // provided in call to this function
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0011", std::to_string(newSeriesNumber));
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0012", AcquisitionNumber);
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0013", InstanceNumber);
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0052", frameOfReferenceUID); // apply
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0018", newFusedSOPInstanceUID);
  itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0080", InstitutionName);

  // get Pixel buffer from fused
  CImageType::Pointer fusedNImage = fused;
  CImageType::PixelContainer *container22 = fusedNImage->GetPixelContainer();
  CImageType::PixelType *buffer22 = container22->GetBufferPointer();

  // now copy all the DICOM tags over
  using DictionaryType = itk::MetaDataDictionary;
  const DictionaryType &dictionaryIn = gdcmImageIO->GetMetaDataDictionary();
  using MetaDataStringType = itk::MetaDataObject<std::string>;
  auto itr = dictionaryIn.Begin();
  auto end = dictionaryIn.End();

  std::string imagePositionPatient; // we might be able to get them this way, but can we set them?
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0032", imagePositionPatient);
  // perhaps we have to use the parsed values to write them again further down?
  double origin3D[3];
  sscanf(imagePositionPatient.c_str(), "%lf\\%lf\\%lf", &(origin3D[0]), &(origin3D[1]), &(origin3D[2]));
  // fprintf(stdout, "image position patient field: %lf, %lf, %lf\n", origin3D[0], origin3D[1], origin3D[2]);

  std::string imageOrientation;
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0037", imageOrientation); // image orientation patient
  double imageOrientationField[6];
  sscanf(imageOrientation.c_str(), "%lf\\%lf\\%lf\\%lf\\%lf\\%lf", &(imageOrientationField[0]), &(imageOrientationField[1]), &(imageOrientationField[2]),
         &(imageOrientationField[3]), &(imageOrientationField[4]), &(imageOrientationField[5]));
  //fprintf(stdout, "reading 0020 0037 as string: %s\n",  imageOrientation.c_str());
  //fprintf(stdout, "parsing 0020 0037 as %f %f %f  %f %f %f\n", imageOrientationField[0], imageOrientationField[1], imageOrientationField[2], imageOrientationField[3], imageOrientationField[4], imageOrientationField[5]);
  // fprintf(stdout, "image orientation field: %lf, %lf, %lf, %lf, %lf, %lf\n", imageOrientationField[0], imageOrientationField[1],
  //        imageOrientationField[2], imageOrientationField[3], imageOrientationField[4], imageOrientationField[5]);

  std::string sliceThicknessString;
  double sliceThickness = 0.0;
  itk::ExposeMetaData<std::string>(dictionarySlice, "0018|0050", sliceThicknessString);
  sscanf(sliceThicknessString.c_str(), "%lf", &sliceThickness);

  std::string imageInstanceString;
  int imageInstance = 0;
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0013", imageInstanceString);
  sscanf(imageInstanceString.c_str(), "%d", &imageInstance);
  // fprintf(stdout, "FOUND INSTANCE: %d\n", imageInstance); // start counting with 0 when we use this value to pick the slice

  std::string sliceLocationString;
  float sliceLocation = 0.0f;
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|1041", sliceLocationString);
  sscanf(sliceLocationString.c_str(), "%f", &sliceLocation);

  std::string imageAcquisitionString;
  int acquisitionNumber = 0;
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0012", imageAcquisitionString);
  sscanf(imageAcquisitionString.c_str(), "%d", &acquisitionNumber);

  // go into slice by applying offset
  // here the problem is that slices can be in a different order from the number in sliceNames
  // we should therefore use an i that corresponds to the image instance number - not the slice name
  // sort order
  uint32_t len = inRegion.GetSize()[0] * inRegion.GetSize()[1] * 3 * sizeof(unsigned char);
  pixeldata.SetByteValue((char *)(((unsigned char *)buffer22)), len);
  image.SetDataElement(pixeldata);

  // create an image (see
  // http://gdcm.sourceforge.net/html/GenFakeImage_8cxx-example.html#_a1)
  gdcm::SmartPointer<gdcm::File> file = new gdcm::File; // empty file
  gdcm::FileDerivation fd;
  const char ReferencedSOPClassUID[] = "1.2.840.10008.5.1.4.1.1.7"; // Secondary Capture
  // create a new frameOfRefenceUID
  gdcm::UIDGenerator fuid;
  fuid.SetRoot("1.3.6.1.4.1.45037");
  // do we need a new one here? - only if we don't find one from before
  if (frameOfReferenceUID == "")
    frameOfReferenceUID = fuid.Generate();

  fd.AddReference(ReferencedSOPClassUID, frameOfReferenceUID.c_str());
  fd.SetPurposeOfReferenceCodeSequenceCodeValue(
      121324); // segmentation  (see
               // https://github.com/malaterre/GDCM/blob/master/Source/MediaStorageAndFileFormat/gdcmFileDerivation.cxx)
  // CID 7203 Image Derivation
  // { "DCM",113072,"Multiplanar reformatting" },
  fd.SetDerivationCodeSequenceCodeValue(113076);
  fd.SetFile(*file);
  // If all Code Value are ok the filter will execute properly
  if (!fd.Derive()) {
    std::cerr << "Sorry could not derive using input info" << std::endl;
    return;
  }
  gdcm::DataSet &ds = fd.GetFile().GetDataSet();
  gdcm::Anonymizer ano;
  ano.SetFile(fd.GetFile());
  std::string seriesDescription;
  int seriesNumber;
  while (itr != end) {
    itk::MetaDataObjectBase::Pointer entry = itr->second;
    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    if (entryvalue) {
      std::string tagkey = itr->first;
      std::string labelId;
      bool found = itk::GDCMImageIO::GetLabelFromTag(tagkey, labelId);
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      if (strcmp(tagkey.c_str(), "0008|103e") == 0) {
        seriesDescription = tagvalue;
      }
      if (strcmp(tagkey.c_str(), "0020|0011") == 0) {
        seriesNumber = atoi(tagvalue.c_str());
      }
      //              if (strcmp(tagkey.c_str(), "0020|1041") == 0) {
      //                // don't overwrite the slice position
      //                ++itr;
      //                continue;
      //              }
      // change window level from -400..600 to 150..180 (why don't we use the computed values? or 0 and 255?)
      if (strcmp(tagkey.c_str(), "0028|1050") == 0) {
        tagvalue = std::string("150");
      }
      if (strcmp(tagkey.c_str(), "0028|1051") == 0) {
        tagvalue = std::string("260");
      }

      unsigned int f1;
      unsigned int f2;
      sscanf(tagkey.c_str(), "%x|%x", &f1, &f2);
      ano.Replace(gdcm::Tag(f1, f2), tagvalue.c_str());
    }
    ++itr;
  }
  gdcm::Attribute<0x0008, 0x2111> at1; // Derivative Description
  at1.SetValue("Fused Segmentation");
  ds.Replace(at1.GetAsDataElement());

  gdcm::Attribute<0x0008, 0x0060> at2; // Derivative Description
  at2.SetValue("MR");
  ds.Replace(at2.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x000E> at3;
  at3.SetValue(newFusedSeriesInstanceUID);
  ds.Replace(at3.GetAsDataElement());

  gdcm::Attribute<0x0008, 0x103E> at4;
  std::string extension = " (fused segmentation)";
  std::ostringstream value;
  value.str("");
  value << seriesDescription;
  // This is a long string and there is a 64 character limit in the
  // standard
  unsigned lengthDesc = value.str().length();
  std::string seriesDesc(value.str(), 0, lengthDesc + extension.length() > 64 ? 64 - extension.length() : lengthDesc + extension.length());
  // itk::EncapsulateMetaData<std::string>(dictionary, "0008|103e", seriesDesc + extension);
  at4.SetValue(seriesDesc + extension);
  ds.Replace(at4.GetAsDataElement());

  // seriesInstance
  gdcm::Attribute<0x0020, 0x0011> at5;
  at5.SetValue(1000 + seriesNumber + 3);
  ds.Replace(at5.GetAsDataElement());

  // use a unique SOPInstanceUID
  gdcm::Attribute<0x0008, 0x0018> at6;
  at6.SetValue(newFusedSOPInstanceUID);
  ds.Replace(at6.GetAsDataElement());

  // set image type to derived
  gdcm::Attribute<0x0008, 0x0008> at_image_type;
  static const gdcm::CSComp values[] = {"DERIVED","SECONDARY"};
  at_image_type.SetValues( values, 2, true ); // true => copy data !
  if ( ds.FindDataElement( at_image_type.GetTag() ) ) {
    const gdcm::DataElement &de = ds.GetDataElement( at_image_type.GetTag() );
    //at_image_type.SetFromDataElement( de );
    // Make sure that value #1 is at least 'DERIVED', so override in all cases:
    at_image_type.SetValue( 0, values[0] );
    at_image_type.SetValue( 1, values[1] );
  }
  ds.Replace( at_image_type.GetAsDataElement() );

  // image position patient from input
  // These values are actually not getting written to the files (RGB has no origin, values are 0\0\0, but see set origin further down)
  gdcm::Attribute<0x0020, 0x0032> at7;
  at7.SetValue(origin3D[0], 0);
  at7.SetValue(origin3D[1], 1);
  at7.SetValue(origin3D[2], 2);
  ds.Replace(at7.GetAsDataElement());
  std::ostringstream value2;
  value2.str("");
  at7.Print(value2);
  // fprintf(stdout, "origin is now supposed to be: %lf\\%lf\\%lf %s\n", origin3D[0], origin3D[1], origin3D[2], value2.str().c_str());
  // For RGB we can set this to make sure they show up at the right location in Horos/OsiriX
  image.SetOrigin(0, origin3D[0]);
  image.SetOrigin(1, origin3D[1]);
  image.SetOrigin(2, origin3D[2]);

  gdcm::Attribute<0x0018, 0x0050> at8;
  at8.SetValue(sliceThickness);
  ds.Replace(at8.GetAsDataElement());

  // TODO: this does not work. We end up with the wrong ImageOrientationPatient information in the output.
  // Either we do not write these attributes for a secondary capture image or we use the correct once 
  // (for linking in Horos for example).
  image.SetDirectionCosines(imageOrientationField);
/*  gdcm::Attribute<0x0020, 0x0037> at9;
  at9.SetValue(imageOrientationField[0], 0);
  at9.SetValue(imageOrientationField[1], 1);
  at9.SetValue(imageOrientationField[2], 2);
  at9.SetValue(imageOrientationField[3], 3);
  at9.SetValue(imageOrientationField[4], 4);
  at9.SetValue(imageOrientationField[5], 5);
  ds.Replace(at9.GetAsDataElement());*/

  // gdcm::Attribute<0x0020, 0x0013> at10;
  // at10.SetValue(imageInstance);
  // ds.Replace(at10.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x1041> at11;
  at11.SetValue(sliceLocation);
  ds.Replace(at11.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x0012> at12;
  at12.SetValue(1000 + acquisitionNumber + 3);
  ds.Replace(at12.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x0013> at13;
  at13.SetValue(imageInstance); // count starts at 1 and increments for all slices
  ds.Replace(at13.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x0052> at14;
  at14.SetValue(frameOfReferenceUID.c_str());
  ds.Replace(at14.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x0010> at16;
  at16.SetValue(studyID.c_str());
  ds.Replace(at16.GetAsDataElement());

  gdcm::Attribute<0x0020, 0x000e> at15;
  at15.SetValue(newFusedSeriesInstanceUID);
  ds.Replace(at15.GetAsDataElement());

  gdcm::ImageWriter writer;
  writer.SetImage(image);
  writer.SetFile(fd.GetFile());
  // std::ostringstream o;
  // o << outputSeries << "/dicom" << i << ".dcm";
  writer.SetFileName(p_out.c_str());
  if (!writer.Write()) {
    return;
  }
  //if (verbose)
  //  fprintf(stdout, "done with writeSecondaryCapture...\n");
  return;
}

bool invalidChar(char c) { return !isprint(static_cast<unsigned char>(c)); }
void stripUnicode(std::string &str) { str.erase(remove_if(str.begin(), str.end(), invalidChar), str.end()); }


template<typename T>
static inline double Lerp(T v0, T v1, T t) {
    return (1 - t)*v0 + t*v1;
}

template<typename T, typename U>
static inline std::vector<T> Quantile(const std::vector<U>& inData, const std::vector<T>& probs) {
    if (inData.empty())
    {
        return std::vector<T>();
    }

    if (1 == inData.size())
    {
        return std::vector<T>(1, inData[0]);
    }

    std::vector<U> data = inData;
    std::sort(data.begin(), data.end());
    std::vector<T> quantiles;

    for (size_t i = 0; i < probs.size(); ++i)
    {
        T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);

        size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
        size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

        T datLeft = data.at(left);
        T datRight = data.at(right);

        T quantile = Lerp<T>(datLeft, datRight, poi - left);

        quantiles.push_back(quantile);
    }

    return quantiles;
}

typedef itk::Image<PixelType, 3> MaskImageType;
typedef itk::Image<PixelType, 3> ImageType3D;

typedef itk::Image<float, 3> InternalImageType;
// typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::Neighborhood<float, 3> NeighborhoodType;
typedef itk::Statistics::DenseFrequencyContainer2 FrequencyContainerType;

typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType3D> Image2CoOccuranceType;
typedef itk::Statistics::ScalarImageToRunLengthMatrixFilter<ImageType3D> Image2RunLengthType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType;
typedef ImageType3D::OffsetType OffsetType;
// typedef itk::AddImageFilter<InternalImageType> AddImageFilterType;
// typedef itk::MultiplyImageFilter<InternalImageType> MultiplyImageFilterType;

using LabelType = unsigned short;
using ShapeLabelObjectType = itk::ShapeLabelObject<LabelType, 3>;

std::map<std::string, std::string> calcTextureRunLengthFeatureImage(ImageType3D::Pointer inputImage, ShapeLabelObjectType *labelObject,
                                                           int pixelMinVal, int pixelMaxVal) {

  ImageType3D::Pointer mask = ImageType3D::New();
  ImageType3D::RegionType maskRegion = inputImage->GetLargestPossibleRegion();
  mask->SetRegions(maskRegion);
  mask->Allocate();
  mask->FillBuffer(itk::NumericTraits<PixelType>::Zero);
  mask->SetOrigin(inputImage->GetOrigin());
  mask->SetSpacing(inputImage->GetSpacing());
  mask->SetDirection(inputImage->GetDirection());
  // use the mask as a single region of interest
  double Np = labelObject->Size();
  for (unsigned int pixelID = 0; pixelID < labelObject->Size(); pixelID++) {
    mask->SetPixel(labelObject->GetIndex(pixelID), 1);
  }

  std::stringstream buf;
  std::map<std::string, std::string> results;

  // Image2RunLengthType
  using FilterType = itk::Statistics::ScalarImageToRunLengthMatrixFilter<ImageType3D>;
  auto filter = FilterType::New();
  filter->SetInput(inputImage);
  ImageType3D::OffsetType      offset1 = { { 0, -1 } };
  ImageType3D::OffsetType      offset2 = { { -1, 0 } };
  FilterType::OffsetVectorPointer offsetV = FilterType::OffsetVector::New();
  offsetV->push_back(offset1);
  offsetV->push_back(offset2);

  filter->SetOffsets(offsetV);
  filter->SetMaskImage(mask);
  // purposely setting the max value to max(Image)+1
  filter->SetPixelValueMinMax(pixelMinVal, pixelMaxVal+1);
  int Nr = 8;
  filter->SetDistanceValueMinMax(0, Nr);
  int Ng = 5;
  filter->SetNumberOfBinsPerAxis(Ng);
  filter->Update();
  const FilterType::HistogramType * hist = filter->GetOutput();
  // resulting frequencies are in hist->GetMeasurementVectorSize() entries
  // from these frequencies we need to compute some summary values
  // see https://pyradiomics.readthedocs.io/en/v1.0/radiomics.html
  int binsPerAxis = filter->GetNumberOfBinsPerAxis();
  using IndexType = FilterType::HistogramType::IndexType;
  double sum = 0.0;
  double sum2 = 0.0;
  for (unsigned int i = 0; i < Ng; i++) {
    for (unsigned int j = 0; j < Nr; j++) {
        IndexType index(hist->GetMeasurementVectorSize());
        index[0] = i;
        index[1] = j;
        buf.str("");
        float val = hist->GetFrequency(index);
        sum += val;
        sum2 += (val*val);
        buf << val;
        char column_name[256];
        snprintf(column_name, 256, "tex_runl_hist_%d_%d", i, j);
        results.insert(std::pair<std::string, std::string>(column_name, buf.str()));
    }
  }

  // in order to compute the features for run length we need the sum (denominator)
  // SRE
  double nominator = 0.0;
  for (unsigned int i = 0; i < Ng; i++) {
    for (unsigned int j = 0; j < Nr; j++) {
        IndexType index(hist->GetMeasurementVectorSize());
        index[0] = i;
        index[1] = j;
        float val = hist->GetFrequency(index);
        nominator += val/((i+1)*(i+1));
    }
  }
  buf.str("");
  buf << nominator / sum;
  results.insert(std::pair<std::string, std::string>("tex_runl_sre", buf.str()));

  // LRE
  nominator = 0.0;
  for (unsigned int i = 0; i < Ng; i++) {
    for (unsigned int j = 0; j < Nr; j++) {
        IndexType index(hist->GetMeasurementVectorSize());
        index[0] = i;
        index[1] = j;
        float val = hist->GetFrequency(index);
        nominator += val*((j+1)*(j+1));
    }
  }
  buf.str("");
  buf << nominator / sum;
  results.insert(std::pair<std::string, std::string>("tex_runl_gre", buf.str()));

  // GLN
  nominator = 0.0;
  for (unsigned int i = 0; i < Ng; i++) {
    double sum_tmp = 0.0;
    for (unsigned int j = 0; j < Nr; j++) {
        IndexType index(hist->GetMeasurementVectorSize());
        index[0] = i;
        index[1] = j;
        float val = hist->GetFrequency(index);
        sum_tmp += val;
    }
    nominator += sum_tmp*sum_tmp;
  }
  buf.str("");
  buf << nominator / sum;
  results.insert(std::pair<std::string, std::string>("tex_runl_gln", buf.str()));

  // GLNN
  buf.str("");
  buf << nominator / sum2;
  results.insert(std::pair<std::string, std::string>("tex_runl_glnn", buf.str()));

  // RLN
  nominator = 0.0;
  for (unsigned int j = 0; j < Nr; j++) {
    double sum_tmp = 0.0;
    for (unsigned int i = 0; i < Ng; i++) {
        IndexType index(hist->GetMeasurementVectorSize());
        index[0] = i;
        index[1] = j;
        float val = hist->GetFrequency(index);
        sum_tmp += val;
    }
    nominator += sum_tmp*sum_tmp;
  }
  buf.str("");
  buf << nominator / sum;
  results.insert(std::pair<std::string, std::string>("tex_runl_rln", buf.str()));

  // RLNN
  buf.str("");
  buf << nominator / sum2;
  results.insert(std::pair<std::string, std::string>("tex_runl_rlnn", buf.str()));

  // RP
  nominator = 0.0;
  for (unsigned int j = 0; j < Nr; j++) {
    for (unsigned int i = 0; i < Ng; i++) {
        IndexType index(hist->GetMeasurementVectorSize());
        index[0] = i;
        index[1] = j;
        float val = hist->GetFrequency(index);
        nominator += val/Np;
    }
  }
  buf.str("");
  buf << nominator;
  results.insert(std::pair<std::string, std::string>("tex_runl_rp", buf.str()));

  return results;
}


std::map<std::string, std::string> calcTextureFeatureImage(OffsetType offset, ImageType3D::Pointer inputImage, ShapeLabelObjectType *labelObject,
                                                           int pixelMinVal, int pixelMaxVal) {
  // Gray Level Co-occurance Matrix Generator
  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
  glcmGenerator->SetOffset(offset);

  ImageType3D::Pointer mask = ImageType3D::New();
  ImageType3D::RegionType maskRegion = inputImage->GetLargestPossibleRegion();
  mask->SetRegions(maskRegion);
  mask->Allocate();
  mask->FillBuffer(itk::NumericTraits<PixelType>::Zero);
  mask->SetOrigin(inputImage->GetOrigin());
  mask->SetSpacing(inputImage->GetSpacing());
  mask->SetDirection(inputImage->GetDirection());
  // use the mask as a single region of interest
  for (unsigned int pixelID = 0; pixelID < labelObject->Size(); pixelID++) {
    mask->SetPixel(labelObject->GetIndex(pixelID), 1);
  }

  glcmGenerator->SetMaskImage(mask);
  // glcmGenerator->SetInsidePixelValue(1);  // default value for inside is 1 - so only if our mask has a single value this will work, if we have more label we
  // need to do this for each label...

  glcmGenerator->SetNumberOfBinsPerAxis(16);                    // reasonable number of bins
  glcmGenerator->SetPixelValueMinMax(pixelMinVal, pixelMaxVal); // for input UCHAR pixel type (but we are using unsigned short with a max value of maybe 311)
  Hist2FeaturesType::Pointer featureCalc = Hist2FeaturesType::New();
  // Region Of Interest
  typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> roiType;
  roiType::Pointer roi = roiType::New();
  roi->SetInput(inputImage);

  /*  InternalImageType::RegionType window;
    InternalImageType::RegionType::SizeType size;
    size.Fill(50);
    window.SetSize(size);

    window.SetIndex(0, 0);
    window.SetIndex(1, 0);
    window.SetIndex(2, 0);*/

  // MaskImageType::RegionType maskRegion = maskImage->GetLargestPossibleRegion();
  roi->SetRegionOfInterest(maskRegion);
  roi->Update();

  glcmGenerator->SetInput(roi->GetOutput());
  glcmGenerator->Update();

  featureCalc->SetInput(glcmGenerator->GetOutput());
  featureCalc->Update();

  std::stringstream buf;
  std::map<std::string, std::string> results;
  buf.str("");
  buf << featureCalc->GetEntropy();
  results.insert(std::pair<std::string, std::string>("tex_entropy", buf.str()));

  buf.str("");
  buf << featureCalc->GetEnergy();
  results.insert(std::pair<std::string, std::string>("tex_energy", buf.str()));

  buf.str("");
  buf << featureCalc->GetCorrelation();
  results.insert(std::pair<std::string, std::string>("tex_correlation", buf.str()));

  buf.str("");
  buf << featureCalc->GetInertia();
  results.insert(std::pair<std::string, std::string>("tex_inertia", buf.str()));

  buf.str("");
  buf << featureCalc->GetHaralickCorrelation();
  results.insert(std::pair<std::string, std::string>("tex_haralick_correlation", buf.str()));

  buf.str("");
  buf << featureCalc->GetInverseDifferenceMoment();
  results.insert(std::pair<std::string, std::string>("tex_inverse_difference_moment", buf.str()));

  buf.str("");
  buf << featureCalc->GetClusterProminence();
  results.insert(std::pair<std::string, std::string>("tex_cluster_prominence", buf.str()));

  buf.str("");
  buf << featureCalc->GetClusterShade();
  results.insert(std::pair<std::string, std::string>("tex_cluster_shade", buf.str()));

  /*  std::cout << "\n Entropy : ";
    std::cout << featureCalc->GetEntropy() << "\n Energy";
    std::cout << featureCalc->GetEnergy() << "\n Correlation";
    std::cout << featureCalc->GetCorrelation() << "\n Inertia";
    std::cout << featureCalc->GetInertia() << "\n HaralickCorrelation";
    std::cout << featureCalc->GetHaralickCorrelation() << "\n InverseDifferenceMoment";
    std::cout << featureCalc->GetInverseDifferenceMoment() << "\nClusterProminence";
    std::cout << featureCalc->GetClusterProminence() << "\nClusterShade";
    std::cout << featureCalc->GetClusterShade(); */

  return results;
}

void computeBiomarkers(Report *report, std::string output_path, std::string imageSeries, std::string labelSeries, bool isMosaic) {
  if (verbose) 
    fprintf(stdout, "start computing biomarkers...\n");
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(false); // we want to use the keys as SeriesInstanceUIDs
  nameGenerator->AddSeriesRestriction("0008|0060");
  nameGenerator->SetRecursive(true);
  nameGenerator->SetDirectory(output_path);

  typedef std::vector<std::string> FileNamesContainer;
  FileNamesContainer fileNames;  // for the label series
  FileNamesContainer fileNames2; // for the image series

  typedef itk::ImageSeriesReader<MaskImageType> MaskReaderType;
  MaskReaderType::Pointer reader = MaskReaderType::New();

  typedef itk::ImageSeriesReader<ImageType3D> ImageReaderType;
  ImageReaderType::Pointer readerImage = ImageReaderType::New();

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  dicomIO->LoadPrivateTagsOn();

  ImageIOType::Pointer dicomIOImage = ImageIOType::New();
  dicomIOImage->LoadPrivateTagsOn();

  reader->SetImageIO(dicomIO);
  readerImage->SetImageIO(dicomIOImage);

  fileNames = nameGenerator->GetFileNames(labelSeries);
  reader->SetFileNames(fileNames);

  try {
    reader->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return;
  }

  fileNames2 = nameGenerator->GetFileNames(imageSeries);
  readerImage->SetFileNames(fileNames2);

  try {
    readerImage->Update();
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return;
  }

  // we should read in the image series as well - we want to compute some biomarkers from that series
  ImageType3D::Pointer image = readerImage->GetOutput();

  using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;

  MaskImageType::Pointer mask = reader->GetOutput();
  // using the above generator we can read in the image and label series again as volumes

  using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<MaskImageType, MaskImageType>;
  using I2LType = itk::LabelImageToShapeLabelMapFilter<MaskImageType, LabelMapType>;

  if (verbose) {
    fprintf(stdout, "compute connected components on mask image...\n"); 
    fflush(stdout);
  }
  auto connected = ConnectedComponentImageFilterType::New();
  connected->SetInput(mask);
  connected->Update();
  if (verbose) {
    fprintf(stdout, "connected components done.\n"); 
    fflush(stdout);
  }

  if (verbose) {
    fprintf(stdout, "compute label to shape filter...\n"); 
    fflush(stdout);
  }
  using I2LType = itk::LabelImageToShapeLabelMapFilter<MaskImageType, LabelMapType>;
  auto i2l = I2LType::New();
  i2l->SetInput(connected->GetOutput());
  i2l->SetComputePerimeter(true);
  i2l->SetComputeFeretDiameter(true);
  i2l->SetComputeOrientedBoundingBox(true);
  i2l->Update();
  if (verbose) {
    fprintf(stdout, "label to shape filter done.\n"); 
    fflush(stdout);
  }

  LabelMapType *labelMap = i2l->GetOutput();

  NeighborhoodType neighborhood;
  neighborhood.SetRadius(1); // 3x3x3
  unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
  OffsetType offset;

  // Retrieve all attributes
  std::stringstream buf;
  int imageMin = 0;
  int imageMax = 0;
  int imageMean = 0;
  if (verbose) {
    fprintf(stdout, "found %ld label%s\n", labelMap->GetNumberOfLabelObjects(), labelMap->GetNumberOfLabelObjects()!= 1?"s":"");
  }

  // compute a key image we can use in the report
  generateImageReturn rets;
  if (isMosaic) {
    rets = generateKeyImageMosaic(image, labelMap, std::vector<int>{512,512}, report->BrightnessContrastLL, report->BrightnessContrastUL);
  } else {
    rets = generateKeyImage(image, labelMap, std::vector<int>{512,512}, report->BrightnessContrastLL, report->BrightnessContrastUL);
  }
  report->keyImage = rets.keyImage;
  report->keyImagePositions = rets.pos;
  report->keyImageTexts = rets.text;

  // need to add some text to this image, number and volume of this location
  // make summary object large enough for all rois
  for (int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n) {
    if (report->summary.size() < n+1)
      report->summary.push_back(std::vector<std::string>{report->summary[0][0]});
    // report->measures
    if (report->measures.size() < n+1) {
      std::map<std::string, std::string> *meas = new std::map<std::string, std::string>();
      report->measures.push_back(*meas);
    }
  }

  // do the summary per label object
  for (unsigned int n_idx = 0; n_idx < labelMap->GetNumberOfLabelObjects(); ++n_idx) {
    unsigned int n = n_idx; // report->roi_order[n_idx];
    //if (report->summary.size() < n_idx+1)
    //  report->summary.push_back(std::vector<std::string>{report->summary[0][0]}); // TODO: add an empty line here, or get the first line of the previous/first summary
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n); // the label number is the connected component number - not the one label as mask
    if (verbose) {
      fprintf(stdout, "processing label %d\n", labelObject->GetLabel()); fflush(stdout);
    }

    std::map<std::string, std::string> *meas = &(report->measures[n_idx]); // new std::map<std::string, std::string>();

    buf.str("");
    buf << itk::NumericTraits<LabelMapType::LabelType>::PrintType(labelObject->GetLabel());
    meas->insert(std::make_pair("region_number", buf.str()));
    buf.str("");
    buf << "3D connected region: " << itk::NumericTraits<LabelMapType::LabelType>::PrintType(labelObject->GetLabel());
    report->summary[n].push_back(buf.str());
    // what image is this based on?
    report->summary[n].push_back("image: " + imageSeries + ", label: " + labelSeries);

    buf.str("");
    buf << labelObject->GetBoundingBox();
    meas->insert(std::make_pair("boundingbox", buf.str()));
    buf.str(""); // clear the buffer
    buf << "    BoundingBox: " << labelObject->GetBoundingBox();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetNumberOfPixels();
    meas->insert(std::make_pair("number_of_pixel", buf.str()));
    buf.str("");
    buf << "    NumberOfPixels: " << labelObject->GetNumberOfPixels();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetPhysicalSize();
    meas->insert(std::make_pair("physical_size", buf.str()));
    buf.str("");
    buf << "    PhysicalSize: " << labelObject->GetPhysicalSize();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetCentroid();
    meas->insert(std::make_pair("centroid", buf.str()));
    buf.str("");
    buf << "    Centroid: " << labelObject->GetCentroid();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetNumberOfPixelsOnBorder();
    meas->insert(std::make_pair("pixel_on_border", buf.str()));
    buf.str("");
    buf << "    NumberOfPixelsOnBorder: " << labelObject->GetNumberOfPixelsOnBorder();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetPerimeterOnBorder();
    meas->insert(std::make_pair("perimeter_on_border", buf.str()));
    buf.str("");
    buf << "    PerimeterOnBorder: " << labelObject->GetPerimeterOnBorder();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetFeretDiameter();
    meas->insert(std::make_pair("feret_diameter", buf.str()));
    buf.str("");
    buf << "    FeretDiameter: " << labelObject->GetFeretDiameter();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetPrincipalMoments();
    meas->insert(std::make_pair("principal_moments", buf.str()));
    buf.str("");
    buf << "    PrincipalMoments: " << labelObject->GetPrincipalMoments();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetPrincipalAxes();
    meas->insert(std::make_pair("principal_axes", buf.str()));
    buf.str("");
    buf << "    PrincipalAxes: " << labelObject->GetPrincipalAxes();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetElongation();
    meas->insert(std::make_pair("elongation", buf.str()));
    buf.str("");
    buf << "    Elongation: " << labelObject->GetElongation();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetPerimeter();
    meas->insert(std::make_pair("perimeter", buf.str()));
    buf.str("");
    buf << "    Perimeter: " << labelObject->GetPerimeter();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetRoundness();
    meas->insert(std::make_pair("roundness", buf.str()));
    buf.str("");
    buf << "    Roundness: " << labelObject->GetRoundness();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetEquivalentSphericalRadius();
    meas->insert(std::make_pair("equivalent_spherical_radius", buf.str()));
    buf.str("");
    buf << "    EquivalentSphericalRadius: " << labelObject->GetEquivalentSphericalRadius();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetEquivalentSphericalPerimeter();
    meas->insert(std::make_pair("equivalent_spherical_perimeter", buf.str()));
    buf.str("");
    buf << "    EquivalentSphericalPerimeter: " << labelObject->GetEquivalentSphericalPerimeter();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetEquivalentEllipsoidDiameter();
    meas->insert(std::make_pair("equivalent_ellipsoid_diameter", buf.str()));
    buf.str("");
    buf << "    EquivalentEllipsoidDiameter: " << labelObject->GetEquivalentEllipsoidDiameter();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetFlatness();
    meas->insert(std::make_pair("flatness", buf.str()));
    buf.str("");
    buf << "    Flatness: " << labelObject->GetFlatness();
    report->summary[n].push_back(buf.str());

    buf.str("");
    buf << labelObject->GetPerimeterOnBorderRatio();
    meas->insert(std::make_pair("perimeter_on_border_ratio", buf.str()));
    buf.str("");
    buf << "    PerimeterOnBorderRatio: " << labelObject->GetPerimeterOnBorderRatio();
    report->summary[n].push_back(buf.str());

    // we should check for this labelObject in image what the intensities are
    if (1) {
      double sum = 0.0f;
      int v;
      itk::Index<3U> index;
      std::vector<int> pixelValues;
      for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
        index = labelObject->GetIndex(pixelId);
        v = image->GetPixel(index);
        if (pixelId == 0) {
          imageMin = v;
          imageMax = v;
        }
        pixelValues.push_back(v);
        sum += v;
        if (v < imageMin)
          imageMin = v;
        if (v > imageMax)
          imageMax = v;
      }
      imageMean = sum / pixelValues.size();
      int median = 0;
      sort(pixelValues.begin(), pixelValues.end());
      int size = pixelValues.size();
      if ((pixelValues.size() % 2) == 0) {
        median = (pixelValues[size / 2 - 1] + pixelValues[size / 2]) / 2;
      } else {
        median = pixelValues[size / 2];
      }
      float stdev = 0.0f;
      double sum2 = 0.0f;
      for (int i = 0; i < pixelValues.size(); i++) {
        sum2 += (pixelValues[i] - imageMean) * (pixelValues[i] - imageMean);
      }
      sum2 /= (pixelValues.size() - 1);
      stdev = std::sqrt(sum2);

      // compute the Q1 and Q3 as well
      auto quartiles = Quantile<double, int>(pixelValues, {0.25, 0.5, 0.75});

      buf.str("");
      buf << "    Intensity in region min: " << imageMin << ", q1: " << quartiles[0] << ", mean: " << imageMean << ", median: " << median
          << ", q3: " << quartiles[2] << ", max: " << imageMax << ", stdev: " << stdev;
      report->summary[n].push_back(buf.str());

      buf.str("");
      buf << quartiles[0];
      meas->insert(std::make_pair("image_intensity_q1", buf.str()));
      // buf.str("");
      // buf << "    Q1 intensity: " << quartiles[0];
      // report->summary.push_back(buf.str());

      buf.str("");
      buf << quartiles[2];
      meas->insert(std::make_pair("image_intensity_q3", buf.str()));
      // buf.str("");
      // buf << "    Q3 intensity: " << quartiles[2];
      // report->summary.push_back(buf.str());

      buf.str("");
      buf << stdev;
      meas->insert(std::make_pair("image_intensity_stdev", buf.str()));
      // buf.str("");
      // buf << "    Stdev intensity: " << stdev;
      // report->summary.push_back(buf.str());

      buf.str("");
      buf << imageMin;
      meas->insert(std::make_pair("image_intensity_min", buf.str()));
      // buf.str("");
      // buf << "    Min intensity: " << imageMin;
      // report->summary.push_back(buf.str());

      buf.str("");
      buf << imageMax;
      meas->insert(std::make_pair("image_intensity_max", buf.str()));
      // buf.str("");
      // buf << "    Max intensity: " << imageMax;
      // report->summary.push_back(buf.str());

      buf.str("");
      buf << imageMean;
      meas->insert(std::make_pair("image_intensity_mean", buf.str()));
      // buf.str("");
      // buf << "    Mean intensity: " << (imageMean);
      // report->summary.push_back(buf.str());

      buf.str("");
      buf << sum;
      meas->insert(std::make_pair("image_intensity_sum", buf.str()));
      buf.str("");
      buf << "    Sum intensity: " << (sum);
      report->summary[n].push_back(buf.str());

      buf.str("");
      buf << median;
      meas->insert(std::make_pair("image_intensity_median", buf.str()));
      // buf.str("");
      // buf << "    Median intensity: " << (median);
      // report->summary.push_back(buf.str());

      // for this labelled object
      for (unsigned int d = 0; d < centerIndex; d++) { // compute this for each d
        offset = neighborhood.GetOffset(d);
        std::map<std::string, std::string> textureFeatures = calcTextureFeatureImage(offset, image, labelObject, imageMin, imageMax);
        // now add these features to the report
        std::stringstream row;
        row.str("");
        row << "    ";
        int counter = 0;
        for (std::map<std::string, std::string>::iterator it = textureFeatures.begin(); it != textureFeatures.end(); ++it) {
          std::stringstream key;
          key << it->first << "_" << std::setfill('0') << std::setw(2) << d;
          meas->insert(std::make_pair(key.str(), it->second));
          row << key.str() << ": " << it->second << " ";
          if ((counter + 1) % 4 == 0) {
            report->summary[n].push_back(row.str());
            row.str("");
            row << "    ";
          }
          counter++;
        }
        if (row.str().size() > 4)
          report->summary[n].push_back(row.str());
      }
      // compute the run length information
      std::map<std::string, std::string> textureFeatures = calcTextureRunLengthFeatureImage(image, labelObject, imageMin, imageMax);
      std::stringstream row;
      row.str("");
      row << "    ";
      int counter = 0;
      for (std::map<std::string, std::string>::iterator it = textureFeatures.begin(); it != textureFeatures.end(); ++it) {
        std::stringstream key;
        key << it->first;
        meas->insert(std::make_pair(key.str(), it->second));
        row << key.str() << ": " << it->second << " ";
        if ((counter + 1) % 4 == 0) {
          report->summary[n].push_back(row.str());
          row.str("");
          row << "    ";
        }
        counter++;
      }
      if (row.str().size() > 4)
        report->summary[n].push_back(row.str());

      //report->measures.push_back(*meas);
    }
  }
  if (verbose) 
    fprintf(stdout, "end computing biomarkers...\n");

}

// maskImage = readMaskImage2D(mask, sliceNr);
typedef itk::Image<unsigned char, 2> MaskSliceImageType;
MaskSliceImageType *readMaskImage2D(std::string filename, int sliceNr) {
  MaskSliceImageType *a;
  return a;
}

static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
  rtrim(s);
  return s;
}
int main(int argc, char *argv[]) {
  setlocale(LC_NUMERIC, "en_US.utf-8");

  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  std::string versionString = std::string("0.0.3.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  command.SetVersion(versionString.c_str());
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("ImageAndMask2Report: Creates a report from an image and mask pair (DICOM format). For all files to be created specify a REPORT_FONT_PATH (export REPORT_FONT_PATH=Monaco.ttf).");
  command.SetCategory("image conversion");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("maskdir", "Directory with input mask DICOM series (0 - background, 1 - foreground).", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for images/, labels/, fused/, and reports/ folder as DICOM. The redcap/ folder contains series folders with output.json files for REDCap imports.", MetaCommand::STRING, true);

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.SetOptionLongTag("SeriesName", "seriesname");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);

  command.SetOption("MaskSeriesName", "m", false, "Select series by series name (if more than one series is present).");
  command.SetOptionLongTag("MaskSeriesName", "maskseriesname");
  command.AddOptionField("MaskSeriesName", "maskseriesname", MetaCommand::STRING, false);

  command.SetOption("Info", "i", false, "Specify an info message that will appear in the report. This option can be used to identify the version/container used for creating the segmentation.");
  command.SetOptionLongTag("Info", "info");
  command.AddOptionField("Info", "info", MetaCommand::STRING, false);

  command.SetOption(
      "UIDFixed", "u", false,
      "If enabled identifiers are stable - will not change for a given input. This allows image series to overwrite each other - assuming that the PACS "
      "supports this overwrite mode. By default the SeriesInstanceUID and SOPInstanceUID values are generated again every time the processing is done.");
  command.SetOptionLongTag("UIDFixed", "uid-fixed");

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  command.SetOption("ReportType", "r", false, "Select report type for key image (mosaic|curvilinear)");
  command.SetOptionLongTag("ReportType", "reporttype");
  command.AddOptionField("ReportType", "value", MetaCommand::STRING, false);

  command.SetOption("BrightnessContrastLL", "d", false, "Set threshold for brightness / contrast based on cummulative histogram lower limit (percentage dark pixel 0.01).");
  command.SetOptionLongTag("BrightnessContrastLL", "brightness-contrast-ll");
  command.AddOptionField("BrightnessContrastLL", "value", MetaCommand::FLOAT, false);

  command.SetOption("BrightnessContrastUL", "b", false, "Set threshold for brightness / contrast based on cummulative histogram upper limit (percentage bright pixel 0.999).");
  command.SetOptionLongTag("BrightnessContrastUL", "brightness-contrast-ul");
  command.AddOptionField("BrightnessContrastUL", "value", MetaCommand::FLOAT, false);

  // float mean_mean = 23.31783, float mean_stds = 4.539313
  command.SetOption("ZScoreMean", "z", false, "Set z-scores mean value (default: 23.31783).");
  command.SetOptionLongTag("ZScoreMean", "z-score-mean");
  command.AddOptionField("ZScoreMean", "value", MetaCommand::FLOAT, false);

  command.SetOption("ZScoreStd", "s", false, "Set z-score standard deviation (default: 4.539313).");
  command.SetOptionLongTag("ZScoreStd", "z-score-std");
  command.AddOptionField("ZScoreStd", "value", MetaCommand::FLOAT, false);

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  verbose = true; //print out everything
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  float mean_mean = 23.31783;
  float mean_stds = 4.539313;
  if (command.GetOptionWasSet("ZScoreMean")) {
    mean_mean = command.GetValueAsFloat("ZScoreMean", "value");
    if (verbose) {
      fprintf(stdout, "read new ZScoreMean from command line: %f\n", mean_mean);
    }
  }
  if (command.GetOptionWasSet("ZScoreStd")) {
    mean_stds = command.GetValueAsFloat("ZScoreStd", "value");
    if (verbose) {
      fprintf(stdout, "read new ZScoreStd from command line: %f\n", mean_stds);
    }
  }
  if (verbose) {
    fprintf(stdout, "USE THESE: %f %f\n", mean_mean, mean_stds);
  }

  float brightness_contrast_ll = 0.01;
  float brightness_contrast_ul = 0.999;
  float brightnesscontrast_ll = brightness_contrast_ll;
  float brightnesscontrast_ul = brightness_contrast_ul;
  if (command.GetOptionWasSet("BrightnessContrastLL")) {
    brightnesscontrast_ll = command.GetValueAsFloat("BrightnessContrastLL", "value");
    if (brightnesscontrast_ll < 0 || brightnesscontrast_ll > 1.0) {
      fprintf(stdout, "Warning: lower brightness values not between 0 and 1. Adjusted to 0.01.\n");
      brightnesscontrast_ll = 0.01;
    }
  }
  if (command.GetOptionWasSet("BrightnessContrastUL")) {
    brightnesscontrast_ul = command.GetValueAsFloat("BrightnessContrastUL", "value");
    if (brightnesscontrast_ul < 0 || brightnesscontrast_ul > 1.0) {
      fprintf(stdout, "Warning: upper brightness values not between 0 and 1. Adjusted to 0.999.\n");
      brightnesscontrast_ul = 0.999;
    }
  }
  if (brightnesscontrast_ul < brightnesscontrast_ll) {
    float tmp = brightnesscontrast_ll;
    brightnesscontrast_ll = brightnesscontrast_ul;
    brightnesscontrast_ul = tmp;
  }
  brightness_contrast_ll = brightnesscontrast_ll;
  brightness_contrast_ul = brightnesscontrast_ul;
  if (verbose) {
    fprintf(stdout, "create report with brightness/contrast settings %.03f %.03f\n", brightness_contrast_ll, brightness_contrast_ul);
  }


  bool isMosaic = false;
  if (command.GetOptionWasSet("ReportType")) {
    std::string bla = command.GetValueAsString("ReportType", "value");
    if (command.GetValueAsString("ReportType", "value") == std::string("mosaic")) {
      if (verbose)
        fprintf(stdout, "Info: Selected mosaic report type\n");
      isMosaic = true;
    } else {
      if (verbose)
        fprintf(stdout, "Info: Selected default report type curvilinear\n");
    }
  }

  bool uidFixedFlag = false;
  if (command.GetOptionWasSet("UIDFixed"))
    uidFixedFlag = true;

  bool seriesIdentifierFlag = false;
  bool maskSeriesIdentifierFlag = false;
  std::string input = command.GetValueAsString("indir");
  std::string mask = command.GetValueAsString("maskdir");
  std::string output = command.GetValueAsString("outdir");

  if (input.size() == 0 || output.size() == 0 || mask.size() == 0) {
    return 1;
  }

  if (command.GetOptionWasSet("SeriesName"))
    seriesIdentifierFlag = true;

  if (command.GetOptionWasSet("MaskSeriesName"))
    maskSeriesIdentifierFlag = true;

  std::string infoMessage = command.GetValueAsString("Info", "info");

  std::string seriesName = command.GetValueAsString("SeriesName", "seriesname");
  std::string maskSeriesName = command.GetValueAsString("MaskSeriesName", "maskseriesname");

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  // to read the PR images we need another ReaderType, they are not images
  // lets call a function that will give us back the polygon information we need to process
  // the image data
  // std::map<std::string, std::string> SOPInstanceUID2SeriesInstanceUID;

  // read in the mask

  //parseForPolygons(input, &storage, &SOPInstanceUID2SeriesInstanceUID, verbose);
  //if (storage.size() == 0) {
  //  fprintf(stderr, "\033[0;31mError\033[0m: No presentation state [PS] files found that contain polylines.\n");
  //  exit(-1);
  //}


  //
  // READ the mask image into an array of slices
  //
  std::vector<MaskImageType2D::Pointer> maskSlices;
  typedef itk::ImageSeriesReader<ImageType> MaskReaderType;
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIOMaskStack = ImageIOType::New();
  dicomIOMaskStack->LoadPrivateTagsOn();
  dicomIOMaskStack->KeepOriginalUIDOn();

  maskReader->SetImageIO(dicomIOMaskStack);

  typedef itk::GDCMSeriesFileNames MaskNamesGeneratorType;
  MaskNamesGeneratorType::Pointer maskNameGenerator = MaskNamesGeneratorType::New();

  maskNameGenerator->SetUseSeriesDetails(false); // we want to use the keys as SeriesInstanceUIDs
  maskNameGenerator->AddSeriesRestriction("0008|0060");
  maskNameGenerator->SetRecursive(true);
  maskNameGenerator->SetDirectory(mask);

  typedef std::vector<std::string> MaskFileNamesContainer;
  MaskFileNamesContainer maskFileNames;
  std::string maskSeriesIdentifier;

  try {
    typedef std::vector<std::string> SeriesIdContainer;

    const SeriesIdContainer &seriesUID = maskNameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

    SeriesIdContainer runTheseMasks;
    if (maskSeriesIdentifierFlag) { // If no optional series identifier
      runTheseMasks.push_back(maskSeriesName);
    } else {
      seriesItr = seriesUID.begin();
      seriesEnd = seriesUID.end();
      while (seriesItr != seriesEnd) {
        runTheseMasks.push_back(seriesItr->c_str());
        ++seriesItr;
      }
    }

    seriesItr = runTheseMasks.begin();
    seriesEnd = runTheseMasks.end();
    while (seriesItr != seriesEnd) {
      maskSeriesIdentifier = seriesItr->c_str();
      ++seriesItr;

      if (verbose) {
        std::cout << "Processing mask series: " << std::endl;
        std::cout << "  " << maskSeriesIdentifier << std::endl;
      }
      fflush(stdout);

      maskFileNames = maskNameGenerator->GetFileNames(maskSeriesIdentifier);

      //  loop over all files in this series
      for (int sliceNr = 0; sliceNr < maskFileNames.size(); sliceNr++) {
        // using ImageType2D = itk::Image<PixelType, 2>;
        typedef itk::ImageFileReader<MaskImageType2D> MaskReader2DType;
        MaskReader2DType::Pointer r = MaskReader2DType::New();
        // typedef itk::GDCMImageIO ImageIOType;
        //  we need to find out what for this image the ReferencedSOPInstanceUID is
        //  only draw the contour on that image
        ImageIOType::Pointer dicomIOMask = ImageIOType::New();
        dicomIOMask->LoadPrivateTagsOn();
        dicomIOMask->KeepOriginalUIDOn();

        r->SetImageIO(dicomIOMask);
        r->SetFileName(maskFileNames[sliceNr]);
        try {
          r->Update();
        } catch (itk::ExceptionObject &err) {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
          return EXIT_FAILURE;
        }
        // now changed the slice we are importing
        maskSlices.push_back(r->GetOutput());
        // ImageType2D::Pointer im2change = r->GetOutput();
      }
    }
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  if (maskSlices.size() == 0) {
    fprintf(stderr, "Error: no mask images found.\n");
    return EXIT_FAILURE;
  }

  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  //typedef itk::Image<unsigned char, 3> MaskImageType;

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO2 = ImageIOType::New();
  dicomIO2->LoadPrivateTagsOn();

  reader->SetImageIO(dicomIO2);

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(false); // we want to use the keys as SeriesInstanceUIDs
  nameGenerator->AddSeriesRestriction("0008|0060");
  nameGenerator->SetRecursive(true);
  nameGenerator->SetDirectory(input);

  try {
    typedef std::vector<std::string> SeriesIdContainer;

    const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    /*while (seriesItr != seriesEnd) {
      std::cout << "Found DICOM Series: ";
      std::cout << std::endl;
      std::cout << "  " << seriesItr->c_str() << std::endl;
      ++seriesItr;
    }*/

    std::string seriesIdentifier;
    std::string PatientName("");
    std::string PatientID("");
    std::string StudyDate("");
    std::string StudyTime("");
    std::string ReferringPhysician("");
    std::string AccessionNumber("");
    std::string StudyID("");

    SeriesIdContainer runThese;
    if (seriesIdentifierFlag) { // If no optional series identifier
      runThese.push_back(seriesName);
    } else {
      seriesItr = seriesUID.begin();
      seriesEnd = seriesUID.end();
      while (seriesItr != seriesEnd) {
        runThese.push_back(seriesItr->c_str());
        ++seriesItr;
      }
    }

    seriesItr = runThese.begin();
    seriesEnd = runThese.end();
    while (seriesItr != seriesEnd) {
      seriesIdentifier = seriesItr->c_str();
      ++seriesItr;
      // do we have this seriesIdentifier in storage? What polygons are part of that?

      if (verbose) {
        std::cout << "Processing series: " << std::endl;
        std::cout << "  " << seriesIdentifier << std::endl;
      }

      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames;

      std::string StudyInstanceUID(""); // keep track of the current StudyInstanceUID

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);
      // we should check now if any of these files SOPInstanceUID appears in our array of polylines (storage)
      // read the images one by one
      if (1) {
        // set as soon as we know what the previous SeriesInstanceUID was
        std::string newSeriesInstanceUID(""); // we can use seriesIdentifier here, set it only once

        // now store the slice as a new series
        std::string newFusedSeriesInstanceUID("");

        // keep the SeriesDescription around so we can use it later
        std::string seriesDescription;

        // remember the last SOPInstanceUID from the list of images
        std::string SOPInstanceUID;

        // remember the last InstitutionName
        std::string InstitutionName;

        std::string StudyDescription;

        // we need to reset some attributes because they are stack specific
        PatientName = "";
        PatientID = "";
        ReferringPhysician = "";
        StudyDate = "";
        StudyTime = "";
        AccessionNumber = "";
        StudyID = "";
        InstitutionName = "";
        StudyDescription = "";

        //  loop over all files in this series
        for (int sliceNr = 0; sliceNr < fileNames.size(); sliceNr++) {
          // using ImageType2D = itk::Image<PixelType, 2>;
          typedef itk::ImageFileReader<ImageType2D> Reader2DType;
          typedef itk::ImageFileWriter<ImageType2D> Writer2DType;
          Reader2DType::Pointer r = Reader2DType::New();
          typedef itk::GDCMImageIO ImageIOType;
          ImageIOType::Pointer dicomIO = ImageIOType::New();
          // ImageIOType::Pointer dicomIOMask = ImageIOType::New();
          dicomIO->LoadPrivateTagsOn();
          dicomIO->KeepOriginalUIDOn();
          // dicomIOMask->LoadPrivateTagsOn();
          // dicomIOMask->KeepOriginalUIDOn();
          //  we need to find out what for this image the ReferencedSOPInstanceUID is
          //  only draw the contour on that image

          r->SetImageIO(dicomIO);
          r->SetFileName(fileNames[sliceNr]);
          try {
            r->Update();
          } catch (itk::ExceptionObject &err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
          }
          // now changed the slice we are importing
          ImageType2D::Pointer im2change = r->GetOutput();

          // get some meta-data from the opened file (find out what polygons are relevant)
          typedef itk::MetaDataDictionary DictionaryType;
          itk::MetaDataDictionary &dictionary = r->GetOutput()->GetMetaDataDictionary();

          // DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();
          std::string SeriesInstanceUID;
          SOPInstanceUID = std::string("");
          std::string seriesNumber;
          itk::ExposeMetaData<std::string>(dictionary, "0020|000e", SeriesInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0008|0018", SOPInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0011", seriesNumber);
          itk::ExposeMetaData<std::string>(dictionary, "0020|000d", StudyInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0008|0080", InstitutionName);
          itk::ExposeMetaData<std::string>(dictionary, "0008|103e", seriesDescription);
          itk::ExposeMetaData<std::string>(dictionary, "0008|1030", StudyDescription);

          // make a copy of this image series in the output/images/ folder
          if (1) {
            // w->SetImageIO(dicomIO); // write the output there
            Writer2DType::Pointer w = Writer2DType::New();
            w->SetInput(im2change);
            // we should have a folder for each image series
            boost::filesystem::path p(fileNames[sliceNr]);
            std::string filename_without_extension = (p.filename().string()).substr(0, (p.filename().string()).find_last_of("."));
            boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "images" + boost::filesystem::path::preferred_separator +
                                            seriesIdentifier + boost::filesystem::path::preferred_separator + filename_without_extension.c_str() + ".dcm";
            if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
              // create the output directory
              create_directories(p_out.parent_path());
            }
            w->SetFileName(p_out.c_str());
            w->SetImageIO(dicomIO);

            try {
              w->Update();
            } catch (itk::ExceptionObject &err) {
              std::cerr << "ExceptionObject caught !" << std::endl;
              std::cerr << err << std::endl;
              return EXIT_FAILURE;
            }
          }

          if (PatientName == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0010|0010", PatientName);
          }
          if (PatientID == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0010|0020", PatientID);
          }
          if (ReferringPhysician == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0008|0090", ReferringPhysician);
          }
          if (StudyDate == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0008|0020", StudyDate);
          }
          if (StudyTime == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0008|0030", StudyTime);
          }
          if (AccessionNumber == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0008|0050", AccessionNumber);
          }
          if (StudyID == "") { // we should use the StudyInstanceUID here, will result in more cases that work back in PACS?
            //itk::ExposeMetaData<std::string>(dictionary, "0020|0010", StudyID);
            StudyID = StudyInstanceUID;
          }
          if (StudyDescription == "") {
            itk::ExposeMetaData<std::string>(dictionary, "0008|1030", StudyDescription);
          }

          if (uidFixedFlag) {
            std::string derivedSeriesInstanceUID(seriesIdentifier);
            std::string endString = ".1";
            if (derivedSeriesInstanceUID.substr(derivedSeriesInstanceUID.size() - 2, 2) == ".1")
              endString = ".2";

            // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
            derivedSeriesInstanceUID = derivedSeriesInstanceUID.substr(0, 64 - 3) + endString;
            newSeriesInstanceUID = derivedSeriesInstanceUID;
          } else {
            if (newSeriesInstanceUID == "") {
              gdcm::UIDGenerator uid;
              uid.SetRoot("1.3.6.1.4.1.45037");
              newSeriesInstanceUID = std::string(uid.Generate());
            } // keep reusing else
          }

          MaskImageType2D::Pointer maskImage = maskSlices[sliceNr];
          //if (polyIds.size() > 0) {
            //maskImage = maskSlices[sliceNr];
            //maskImage = readMaskImage2D(mask, sliceNr);
            // once we have masks we want to save a fused image series (masks ontop of im2change)
          //} else { // without polygon just return an empty image
            // fill the input image with 1 (label)
          //  im2change->FillBuffer(itk::NumericTraits<PixelType>::Zero);
          //}

          // InputImageType::Pointer inputImage = reader->GetOutput();
          ImageType2D::RegionType region;
          region = im2change->GetBufferedRegion();
          ImageType2D::SizeType size = region.GetSize();
          // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

          ImageType2D::PixelContainer *container;
          container = im2change->GetPixelContainer();
          container->SetContainerManageMemory(false);
          unsigned int bla = sizeof(ImageType2D::PixelType);
          ImageType2D::PixelType *buffer2 = container->GetBufferPointer();

          ImageType2D::Pointer nImage;
          MaskImageType2D::PixelContainer *container2 = maskImage->GetPixelContainer();
          MaskImageType2D::PixelType *buffer3 = container2->GetBufferPointer();

          // Here we copy all values over, that is 0, 1, 2, 3 but also additional labels
          // that have been selected before (air in intestines for example).
          // buffer2 is destination, buffer3 is source (mask)
          memcpy(buffer2, &(buffer3[0]), size[0] * size[1] * bla);
          // We can clean the data (remove all other label).
          /*for (int k = 0; k < size[0] * size[1]; k++) {
            if (buffer2[k] > 3) {
              buffer2[k] = 0; // set to background
            }
          }*/

          // dilate and erode the im2change (mask) (would be better if we do this in 3D)
          using StructuringElementType = itk::BinaryBallStructuringElement<PixelType, 2>;
          using ErodeFilterType = itk::BinaryErodeImageFilter<ImageType2D, ImageType2D, StructuringElementType>;
          using DilateFilterType = itk::BinaryDilateImageFilter<ImageType2D, ImageType2D, StructuringElementType>;
          DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
          ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

          StructuringElementType structuringElement;
          structuringElement.SetRadius(1); // 3x3 structuring element
          structuringElement.CreateStructuringElement();
          binaryDilate->SetKernel(structuringElement);
          binaryErode->SetKernel(structuringElement);
          binaryErode->SetForegroundValue(1);
          binaryErode->SetBackgroundValue(0);
          binaryDilate->SetForegroundValue(1);
          binaryDilate->SetBackgroundValue(0);
          binaryDilate->SetInput(im2change); // should this be maskImage?
          binaryErode->SetInput(binaryDilate->GetOutput());
          binaryDilate->SetDilateValue(1);
          binaryErode->SetErodeValue(1);

          try {
            binaryErode->Update();
          } catch (itk::ExceptionObject &err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
          }
          // per slice we should compute volumes here
          float intersliceThickness = 0;

          // create a fused image using the mask in binaryErode->GetOutput()
          if (1) {
            if (uidFixedFlag) {
              std::string derivedFusedSeriesInstanceUID(seriesIdentifier);
              std::string endString = ".3";
              if (derivedFusedSeriesInstanceUID.substr(derivedFusedSeriesInstanceUID.size() - 2, 2) == ".3")
                endString = ".4";

              // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
              derivedFusedSeriesInstanceUID = derivedFusedSeriesInstanceUID.substr(0, 64 - 3) + endString;
              newFusedSeriesInstanceUID = derivedFusedSeriesInstanceUID;
            } else { // we can only do this once!!! not in a loop for each slice
              if (newFusedSeriesInstanceUID == "") {
                gdcm::UIDGenerator uid;
                uid.SetRoot("1.3.6.1.4.1.45037");
                newFusedSeriesInstanceUID = std::string(uid.Generate());
              }
            }

            std::string newFusedSOPInstanceUID("");
            if (uidFixedFlag) {
              std::string derivedFusedSOPInstanceUID(SOPInstanceUID);
              std::string endString = ".4";
              if (derivedFusedSOPInstanceUID.substr(derivedFusedSOPInstanceUID.size() - 2, 2) == ".4")
                endString = ".5";

              // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
              derivedFusedSOPInstanceUID = derivedFusedSOPInstanceUID.substr(0, 64 - 3) + endString;
              newFusedSOPInstanceUID = derivedFusedSOPInstanceUID;
            } else {
              gdcm::UIDGenerator uid;
              uid.SetRoot("1.3.6.1.4.1.45037");
              newFusedSOPInstanceUID = std::string(uid.Generate());
            }

            boost::filesystem::path p(fileNames[sliceNr]);
            std::string filename_without_extension = (p.filename().string()).substr(0, (p.filename().string()).find_last_of("."));
            boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "fused" + boost::filesystem::path::preferred_separator +
                                            newFusedSeriesInstanceUID + boost::filesystem::path::preferred_separator + filename_without_extension.c_str() + ".dcm";
            if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
              // create the output directory
              create_directories(p_out.parent_path());
            }
            // it would be cool to add a filtered version of the labels as well, but that works only for a single label...  or?
            // We don't have a split into different regions of interest, we should use connected components here as well
            // we do that already in computeBiomarkers.

            // in case we need to split the regions into connected components do this:
            // but this does not make sense in 2D, our labels must be computed in 3D.
            /*using ConnectedComponents = itk::ConnectedComponentImageFilter<ImageType2D, ImageType2D>;
            ConnectedComponents::Pointer con = ConnectedComponents::New();
            con->SetInput(binaryErode->GetOutput());
            con->SetBackgroundValue(0);
            con->Update(); */
            writeSecondaryCapture(binaryErode->GetOutput(), fileNames[sliceNr], std::string(p_out.c_str()), uidFixedFlag, newFusedSeriesInstanceUID,
                                  newFusedSOPInstanceUID, verbose, brightness_contrast_ll, brightness_contrast_ul);
            // here is a good place to extract some measures from the masked image (mean, min, max, median, sum, intensity, histogram?)
          }

          // copy the values back to the im2change buffer
          ImageType2D::Pointer cleanMask = binaryErode->GetOutput();
          ImageType2D::PixelContainer *container3;
          container3 = cleanMask->GetPixelContainer();
          ImageType2D::PixelType *buffer4 = container3->GetBufferPointer();
          memcpy(buffer2, &(buffer4[0]), size[0] * size[1] * bla);

          // now change something to make a new copy of that file
          int newSeriesNumber = 1000 + atoi(seriesNumber.c_str()) + 1;
          std::string newSOPInstanceUID = std::string("");
          if (uidFixedFlag) {
            newSOPInstanceUID = SOPInstanceUID;
            // fprintf(stderr, "old SOPInstanceUID: %s\n", newSOPInstanceUID.c_str());
            //  end
            std::string endString = ".1";
            if (newSOPInstanceUID.substr(newSOPInstanceUID.size() - 2, 2) == ".1")
              endString = ".2";
            newSOPInstanceUID = newSOPInstanceUID.substr(0, 64 - 3) + endString;
            // fprintf(stderr, "new SOPInstanceUID: %s\n", newSOPInstanceUID.c_str());
          } else {
            gdcm::UIDGenerator uid;
            uid.SetRoot("1.3.6.1.4.1.45037");
            newSOPInstanceUID = uid.Generate();
          }

          dicomIO->KeepOriginalUIDOn();
          itk::MetaDataDictionary &dictionarySlice = r->GetOutput()->GetMetaDataDictionary();
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0011", std::to_string(newSeriesNumber));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0018", newSOPInstanceUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000e", newSeriesInstanceUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0010", StudyInstanceUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0080", InstitutionName);

          // set to DERIVED\\SECONDARY
/*          gdcm::File &fd = r->GetFile();
          gdcm::DataSet &ds = fd.GetDataSet();

          gdcm::Attribute<0x0008,0x0008> at_image_type;
          static const gdcm::CSComp values[] = {"DERIVED","SECONDARY"};
          at_image_type.SetValues( values, 2, true ); // true => copy data !
          if ( ds.FindDataElement( at_image_type.GetTag() ) ) {
            const gdcm::DataElement &de = ds.GetDataElement( at_image_type.GetTag() );
            at_image_type.SetFromDataElement( de );
            // Make sure that value #1 is at least 'DERIVED', so override in all cases:
            at_image_type.SetValue( 0, values[0] );
          }
          ds.Replace( at_image_type.GetAsDataElement() );
*/
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0008", std::string("DERIVED\\SECONDARY"));

          // set the series description (max 64 characters)
          if (seriesDescription != "")
            seriesDescription += " ";
          std::string newSeriesDescription = seriesDescription + "(mask)";
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|103e", newSeriesDescription.substr(0, 64));
          // set window center and window width
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0028|1050", std::to_string(0.5));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0028|1051", std::to_string(1));

          Writer2DType::Pointer w = Writer2DType::New();
          w->SetInput(im2change);
          // create the output filename
          // we should have a folder for each image series
          boost::filesystem::path p(fileNames[sliceNr]);
          std::string filename_without_extension = (p.filename().string()).substr(0, (p.filename().string()).find_last_of("."));
          boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "labels" + boost::filesystem::path::preferred_separator +
                                          newSeriesInstanceUID.c_str() + boost::filesystem::path::preferred_separator + filename_without_extension.c_str() + ".dcm";
          if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
            // fprintf(stderr, "create directory with name: \"%s\" for newSeriesInstanceUID: \"%s\"\n", p_out.c_str(), newSeriesInstanceUID.c_str());
            create_directories(p_out.parent_path());
          }
          w->SetFileName(p_out.c_str());
          w->SetImageIO(dicomIO);

          if (verbose) {
            fprintf(stdout, "write: %s with %s %s\n", p_out.c_str(), newSOPInstanceUID.c_str(), newSeriesInstanceUID.c_str());
          }

          try {
            w->Update();
          } catch (itk::ExceptionObject &err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
          }
          // example: ./Modules/Filtering/ImageIntensity/test/itkPolylineMask2DImageFilterTest.cxx
        }
        // we should remember this mapping in a csv file
        boost::filesystem::path csv_out = output + boost::filesystem::path::preferred_separator + "data.csv";
        if (!itksys::SystemTools::FileIsDirectory(csv_out.parent_path().c_str())) {
          create_directories(csv_out.parent_path());
        }
        if (!boost::filesystem::exists(csv_out)) {
          FILE *fp = fopen(csv_out.c_str(), "a");
          fprintf(fp, "ImageSeriesInstanceUID,LabelSeriesInstanceUID\n");
          fclose(fp);
        }
        FILE *fp = fopen(csv_out.c_str(), "a");
        fprintf(fp, "\"images/%s\",\"labels/%s\"\n", seriesIdentifier.c_str(), newSeriesInstanceUID.c_str());
        fclose(fp);


        // create a report for this series as well
        Report *report = getDefaultReportStruct();
        report->summary = {{"Research PACS Report " + to_simple_string(timeLocal)}};
        // report->summary.push_back(resultJSON["wall_time"]);
        report->StudyInstanceUID = StudyInstanceUID; // but use a new series and SOPInstanceUID generated by getDefaultReportStruct
        report->VersionString = infoMessage; // might be specified by the user or empty string
        report->PatientName = PatientName;
        report->PatientID = PatientID;
        report->StudyDate = StudyDate;
        report->StudyTime = StudyTime;
        report->AccessionNumber = AccessionNumber;
        report->StudyDescription = StudyDescription;
        report->StudyID = StudyID;
        report->SeriesDescription = seriesDescription + " (report)";
        report->ReferringPhysician = ReferringPhysician;
        report->ReportSeriesInstanceUID = seriesIdentifier;
        report->SOPInstanceUID = SOPInstanceUID;
        report->InstitutionName = InstitutionName;
        report->BrightnessContrastLL = brightness_contrast_ll;
        report->BrightnessContrastUL = brightness_contrast_ul;
        report->ReportType = (isMosaic?"mosaic":"curvilinear");

        // TODO: in case we do uid-fixed we would need to create the same report SOPInstanceUID and SeriesInstanceUID
        if (uidFixedFlag) {
          std::string derivedSeriesInstanceUID(seriesIdentifier);
          std::string endString = ".7";
          if (derivedSeriesInstanceUID.substr(derivedSeriesInstanceUID.size() - 2, 2) == ".7")
            endString = ".8";
          // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
          derivedSeriesInstanceUID = derivedSeriesInstanceUID.substr(0, 64 - 3) + endString;
          report->SeriesInstanceUID = derivedSeriesInstanceUID;
          report->ReportSeriesInstanceUID = derivedSeriesInstanceUID;

          std::string newSOPInstanceUID = SOPInstanceUID;
          if (newSOPInstanceUID.substr(newSOPInstanceUID.size() - 2, 2) == ".7")
            endString = ".8";
          report->SOPInstanceUID = newSOPInstanceUID.substr(0, 64 - 3) + endString;
        }
        fflush(stdout);
        // we got the label as a mask stored in the labels folder, read, convert to label and create summary statistics
        computeBiomarkers(report, output, seriesIdentifier, newSeriesInstanceUID, isMosaic);

        int key_fact = 0;
        for (int i = 0; i < report->measures.size(); i++) {
          key_fact += std::stof(report->measures[i].find("physical_size")->second);
        }
        report->key_fact = std::to_string(key_fact);
        report->key_unit = std::string("mm^3");

        // overwrite some report values
        boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "reports" + boost::filesystem::path::preferred_separator +
                                        newSeriesInstanceUID.c_str();
        if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
          if (verbose)
            fprintf(stderr, "create directory with name: \"%s\" for newSeriesInstanceUID: \"%s\"\n", p_out.c_str(), newSeriesInstanceUID.c_str());
          create_directories(p_out.parent_path());
        }
        report->filename = std::string(p_out.c_str());

        // mark a report as empty in case there are no detected regions of interest
        if (!report->keyImage || report->measures.size() == 0) {
          // add to the summary
          if (report->summary.size() > 0) {
            report->summary[0].push_back(std::string(""));
            report->summary[0].push_back(std::string("Empty report, no region of interest could be detected."));
          }
        }

        // default values only make sense if used for the spine segmentation project
        saveReport(report, mean_mean, mean_stds, verbose);

        // add measures to json output, make sure to keep values from previous iteration
        if (!resultJSON.contains("measures")) {
          resultJSON["measures"] = json::array();
        }
        std::vector< std::map<std::string, std::string> >::iterator iter = report->measures.begin();
        for(iter; iter < report->measures.end(); iter++) {
          resultJSON["measures"].push_back(*iter);
        }
        if (!resultJSON.contains("meta-data")) {
          resultJSON["meta-data"] = json::array();
        }
        json obj = { 
          { "PatientName", report->PatientName },
          { "PatientID", report->PatientID },
          { "SeriesInstanceUID", report->SeriesInstanceUID },
          { "StudyInstanceUID", report->StudyInstanceUID },
          { "StudyDescription", report->StudyDescription },
          { "FrameOfReferenceUID", report->FrameOfReferenceUID },
          { "ReferringPhysician", report->ReferringPhysician },
          { "StudyID", report->StudyID },
          { "AccessionNumber", report->AccessionNumber },
          { "SeriesDescription", report->SeriesDescription },
          { "StudyDate", report->StudyDate },
          { "StudyTime", report->StudyTime },
          { "InputImageSeriesIdentifier", seriesIdentifier },
          { "InputMaskSeriesIdentifier", maskSeriesIdentifier },
      	  { "ReportSOPInstanceUID", report->SOPInstanceUID },
          { "InstitutionName", report->InstitutionName },
          { "Summary", report->summary }
        };
        resultJSON["meta-data"].push_back(obj);
        // produce a REDCap friendly output format for the measures
        json redcap = json::array();
        for (int i = 0; i < report->measures.size(); i++) {
          auto m = report->measures[i];
          // go through all the entries in that map
          for (std::map<std::string, std::string>::iterator iter = m.begin(); iter != m.end(); ++iter) {
            std::string key = iter->first;
            std::string value = iter->second;
            json a = json::object();
            a["record_id"] = rtrim_copy(PatientID);
            std::string ref = ReferringPhysician;
            std::string startString("EventName:");
            if (ref.find(startString) == 0) {
              ref = ref.substr(startString.length());
            }
            a["redcap_event_name"] = ref;
            a["value"] = value;
            a["field_name"] = key;
            a["redcap_repeat_instrument"] = "pr2mask";
            a["redcap_repeat_instance"] = std::to_string(i + 1); // start counting with 1 for redcap
            redcap.push_back(a);
          }
        }

        boost::filesystem::path output_out = output + boost::filesystem::path::preferred_separator + "redcap" + boost::filesystem::path::preferred_separator +
                                             newSeriesInstanceUID.c_str() + boost::filesystem::path::preferred_separator + "output.json";
        if (!itksys::SystemTools::FileIsDirectory(output_out.parent_path().c_str())) {
          create_directories(output_out.parent_path());
        }
        std::ofstream out2(output_out.c_str());
        std::string res2 = redcap.dump(4) + "\n";
        out2 << res2;
        out2.close();

        // computational time
        boost::posix_time::ptime timeLocalEnd = boost::posix_time::microsec_clock::local_time();
        boost::posix_time::time_period tp(timeLocal, timeLocalEnd);
        resultJSON["wall_time"] = boost::posix_time::to_simple_string(timeLocalEnd - timeLocal);
        std::string res = resultJSON.dump(4) + "\n";
        // save the json information to a file as well, use folder names
        boost::filesystem::path json_out = output + boost::filesystem::path::preferred_separator + seriesIdentifier + "_" + newSeriesInstanceUID + ".json";
        std::ofstream out(json_out.c_str());
        out << res;
        out.close();

        // for each output.json we should also export a data dictionary on how to interprete this file
        // best would be a zip file that looks like the one exported from REDCap
        // https://gist.github.com/clalancette/bb5069a09c609e2d33c9858fcc6e170e
        // create OriginID.txt and instrument.csv, zip the content of that folder
        boost::filesystem::path originid_file = output + boost::filesystem::path::preferred_separator + "OriginID.txt";
        boost::filesystem::path instrument_file = output + boost::filesystem::path::preferred_separator + "instrument.csv";
        std::ofstream out3(originid_file.c_str());
        std::string content("PR2MASK");
        out3 << content;
        out3.close();

        std::ofstream out4(instrument_file.c_str());
        content = std::string("\"Variable / Field Name\",\"Form Name\",\"Section Header\",\"Field Type\",\"Field Label\",\"Choices, Calculations, OR Slider Labels\",\"Field Note\",\"Text Validation Type OR Show Slider Number\",\"Text Validation Min\",\"Text Validation Max\",Identifier?,\"Branching Logic (Show field only if...)\",\"Required Field?\",\"Custom Alignment\",\"Question Number (surveys only)\",\"Matrix Group Name\",\"Matrix Ranking?\",\"Field Annotation\"\n");
        out4 << content;
        // physical_size,pr2mask,,text,"Physical size",,,number,,,,,,,,,,
        for (int meas = 0; meas < report->measures.size(); meas++) {
          auto m = report->measures[meas];
          // go through all the entries in that map
          for (std::map<std::string, std::string>::iterator iter = m.begin(); iter != m.end(); ++iter) {
            std::string key = iter->first;
            std::string type("");
            content = key + std::string(",pr2mask,,text,\"") + key + "\",,," + type + ",,,,,,,,,,\n";
            out4 << content;
          }
          break;
        }
        out4.close();

        // we should add those two files into a zip file (and delete them again)
        int errorp;
        boost::filesystem::path dictionary_file = output + boost::filesystem::path::preferred_separator + "redcap" + boost::filesystem::path::preferred_separator +
                                             newSeriesInstanceUID.c_str() + boost::filesystem::path::preferred_separator + "output_data_dictionary.zip";
        // if the zip file already exists delete it first
        if (boost::filesystem::exists(dictionary_file)) {
          boost::filesystem::remove(dictionary_file);
        }

        zip_t *zipper = zip_open(dictionary_file.string().c_str(), ZIP_CREATE | ZIP_EXCL, &errorp);
        if (zipper == nullptr) {
          zip_error_t ziperror;
          zip_error_init_with_code(&ziperror, errorp);
          throw std::runtime_error(std::string("Failed to open output file ") + dictionary_file.string() + std::string(": ") + zip_error_strerror(&ziperror));
        }
        // add first file
        zip_source_t *source = zip_source_file(zipper, originid_file.c_str(), 0, 0);
        if (source == nullptr) {
          throw std::runtime_error("Failed to add file to zip: " + std::string(zip_strerror(zipper)));
        }
        if (zip_file_add(zipper, "OriginID.txt", source, ZIP_FL_ENC_UTF_8) < 0) {
          zip_source_free(source);
          throw std::runtime_error("Failed to add file to zip: " + std::string(zip_strerror(zipper)));
        }

        // add second file
        source = zip_source_file(zipper, instrument_file.c_str(), 0, 0);
        if (source == nullptr) {
          throw std::runtime_error("Failed to add file to zip: " + std::string(zip_strerror(zipper)));
        }
        if (zip_file_add(zipper, "instrument.csv", source, ZIP_FL_ENC_UTF_8) < 0) {
          zip_source_free(source);
          throw std::runtime_error("Failed to add file to zip: " + std::string(zip_strerror(zipper)));
        }

        zip_close(zipper);
        if (boost::filesystem::exists(originid_file)) {
          boost::filesystem::remove(originid_file);
        }
        if (boost::filesystem::exists(instrument_file)) {
          boost::filesystem::remove(instrument_file);
        }

      }

    } // loop over series

  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  if (verbose) {
    std::string res = resultJSON.dump(4) + "\n";
    fprintf(stdout, "%s", res.c_str());
  }

  return EXIT_SUCCESS;
}
