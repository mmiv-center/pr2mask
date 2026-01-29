/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkPolylineMask2DScanlineImageFilter_hxx
#define itkPolylineMask2DScanlineImageFilter_hxx

#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkLineIterator.h"
#include "itkPathIterator.h"
#include "itkPolylineMask2DScanlineImageFilter.h"
#include "itkProgressReporter.h"
#include "itkVectorContainer.h"

namespace itk {
/**
 * Constructor
 */
template <typename TInputImage, typename TPolyline, typename TOutputImage>
PolylineMask2DScanlineImageFilter<TInputImage, TPolyline, TOutputImage>::PolylineMask2DScanlineImageFilter() {
  this->SetNumberOfRequiredInputs(2);
}

/**
 *
 */
template <typename TInputImage, typename TPolyline, typename TOutputImage>
void PolylineMask2DScanlineImageFilter<TInputImage, TPolyline, TOutputImage>::SetInput1(const TInputImage *input)

{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, const_cast<TInputImage *>(input));
}

/**
 *
 */
template <typename TInputImage, typename TPolyline, typename TOutputImage>
void PolylineMask2DScanlineImageFilter<TInputImage, TPolyline, TOutputImage>::SetInput2(const TPolyline *input) {
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(1, const_cast<TPolyline *>(input));
}

struct _Point2D {
  float x;
  float y;
};

// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(_Point2D p, _Point2D q, _Point2D r) {
  if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) && q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
    return true;

  return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(_Point2D p, _Point2D q, _Point2D r) {
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

  if (std::abs(val) < 1e-6)
    return 0; // collinear

  return (val > 0) ? 1 : 2; // clock or counter clockwise
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(_Point2D p1, _Point2D q1, _Point2D p2, _Point2D q2) {
  // Find the four orientations needed for general and
  // special cases
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);

  // General case
  if (o1 != o2 && o3 != o4)
    return true;

  // Special Cases
  // p1, q1 and p2 are collinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1, p2, q1))
    return true;

  // p1, q1 and q2 are collinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1, q2, q1))
    return true;

  // p2, q2 and p1 are collinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2, p1, q2))
    return true;

  // p2, q2 and q1 are collinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2, q1, q2))
    return true;

  return false; // Doesn't fall in any of the above cases
}

/**
 *
 */
template <typename TInputImage, typename TPolyline, typename TOutputImage>
void PolylineMask2DScanlineImageFilter<TInputImage, TPolyline, TOutputImage>::GenerateData() {
  using LineIteratorType = LineIterator<TOutputImage>;
  using ImageLineIteratorType = ImageLinearIteratorWithIndex<TOutputImage>;

  using InputImageConstIteratorType = ImageRegionConstIterator<TInputImage>;

  using ImageIndexType = typename TOutputImage::IndexType;
  using PixelType = typename TOutputImage::PixelType;
  using OutputImageIteratorType = ImageRegionIterator<TOutputImage>;

  using VertexType = typename TPolyline::VertexType;
  using VertexListType = typename TPolyline::VertexListType;

  typename TInputImage::ConstPointer inputImagePtr(dynamic_cast<const TInputImage *>(this->ProcessObject::GetInput(0)));
  typename TPolyline::ConstPointer polylinePtr(dynamic_cast<const TPolyline *>(this->ProcessObject::GetInput(1)));
  typename TOutputImage::Pointer outputImagePtr(dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(0)));

  outputImagePtr->SetOrigin(inputImagePtr->GetOrigin());
  outputImagePtr->SetSpacing(inputImagePtr->GetSpacing());
  outputImagePtr->SetDirection(inputImagePtr->GetDirection());
  outputImagePtr->SetRequestedRegion(inputImagePtr->GetRequestedRegion());
  outputImagePtr->SetBufferedRegion(inputImagePtr->GetBufferedRegion());
  outputImagePtr->SetLargestPossibleRegion(inputImagePtr->GetLargestPossibleRegion());
  outputImagePtr->Allocate();
  outputImagePtr->FillBuffer(itk::NumericTraits<PixelType>::Zero);

  PixelType zero_val = NumericTraits<PixelType>::ZeroValue();
  auto u_val = static_cast<PixelType>(0);
  auto b_val = static_cast<PixelType>(2);
  auto f_val = static_cast<PixelType>(255);
  float halfAPixelDirX = inputImagePtr->GetSpacing()[0] / 2.0f;
  float halfAPixelDirY = inputImagePtr->GetSpacing()[1] / 2.0f;

  const VertexListType *container = polylinePtr->GetVertexList();

  typename VertexListType::ConstIterator piter = container->Begin();

  // use a scanline method to set masks

  VertexType startVertex;
  VertexType endVertex;
  VertexType pstartVertex;
  VertexType tmpVertex;
  VertexType tmpVertexA;
  VertexType tmpVertexB;

  // bounding box in pixel
  auto size = outputImagePtr->GetLargestPossibleRegion().GetSize();
  int startx = size[0];
  int starty = size[1];
  int endx = 0;
  int endy = 0;

  /* Check if the polyline coordinates are within the input image and get the bounding box for scanline */
  while (piter != container->End()) {
    tmpVertex = piter.Value();
    const auto tmpIndex = outputImagePtr->TransformPhysicalPointToIndex(tmpVertex);
    if (!outputImagePtr->GetBufferedRegion().IsInside(tmpIndex)) {
      itkExceptionMacro(<< "Polyline vertex is out of bounds (Vertex,Index): " << tmpVertex << ", " << tmpIndex);
    }
    if (startx > tmpIndex[0])
      startx = tmpIndex[0];
    if (starty > tmpIndex[1])
      starty = tmpIndex[1];

    if (endx < tmpIndex[0])
      endx = tmpIndex[0];
    if (endy < tmpIndex[1])
      endy = tmpIndex[1];

    ++piter;
  }
  // define a region for this polygon
  using ImageSizeType = typename TOutputImage::SizeType;
  ImageSizeType regionSize;

  // make the mask one bigger so we capture the beginning
  startx--;
  starty--;
  endx++;
  endy++;
  if (endx > size[0] - 1)
    endx = size[0] - 1;
  if (endy > size[1] - 1)
    endy = size[1] - 1;
  if (startx < 0)
    startx = 0;
  if (starty < 0)
    starty = 0;

  // what about rounding? we compute here with int ... but shouldn't we do a float rounding down and up?
  regionSize[0] = endx - startx; // size along X
  regionSize[1] = endy - starty; // size along Y

  using ImageIndexType = typename TOutputImage::IndexType;
  ImageIndexType start;
  start[0] = startx;
  start[1] = starty;

  using ImageRegionType = typename TOutputImage::RegionType;
  ImageRegionType region;
  region.SetIndex(start);
  region.SetSize(regionSize);
  // fprintf(stderr, "start: %d %d\n", start[0], start[1]);
  // fprintf(stderr, "size: %d %d\n", regionSize[0], regionSize[1]);

  ImageLineIteratorType imitLine(outputImagePtr, region);
  imitLine.SetDirection(0);

  // we need to go through the region
  imitLine.GoToBegin();
  while (!imitLine.IsAtEnd()) { // it matters where we place our lines, inside a pixel or between pixels
                                // do not hit again if you are in the same row -  hit already list?
    // fprintf(stderr, "anfang...\n");
    //   walk the line
    auto indexA = imitLine.GetIndex(); // we should not use the index here but the float coordinates, center at middle of pixel
    using PT = typename TOutputImage::PointType;
    PT floatIndexA;
    outputImagePtr->TransformIndexToPhysicalPoint(indexA, floatIndexA); // hope that is in the middle of the pixel now
    floatIndexA[0] -= halfAPixelDirX;                                   // half a pixel in x direction
    //floatIndexA[1] -= halfAPixelDirY / 4;                               // prevent intersections at end points
    floatIndexA[1] -= halfAPixelDirY;                               // prevent intersections at end points
    // fprintf(stderr, " POINT TO FLOAT: %ld %ld %f %f\n", indexA[0], indexA[1], floatIndexA[0], floatIndexA[1]);
    //  auto indexA2 = outputImagePtr->TransformPhysicalPointToContinuousIndex(tmpVertex);
    ++imitLine;
    bool outside = true;            // every time we start a new line we are outside
    std::vector<int> lineHitMemory; // for every line we can hit a vertex only once, need an iterator that is doing consistent indexing
    lineHitMemory.clear();
    // fprintf(stderr, "start a line...\n");
    while (!imitLine.IsAtEndOfLine()) {
      // what is the location?
      auto indexB = imitLine.GetIndex();
      PT floatIndexB;
      outputImagePtr->TransformIndexToPhysicalPoint(indexB, floatIndexB); // hope that is in the middle of the pixel now
      floatIndexB[0] -= halfAPixelDirX;
      // floatIndexB[1] -= halfAPixelDirY / 4;
      floatIndexB[1] -= halfAPixelDirY;
      // we have now a tiny line between indexA and indexB
      // now we can check the intersection with all the lines we have in the polygon
      tmpVertexA = container->ElementAt(0);
      auto indexLineA = outputImagePtr->TransformPhysicalPointToIndex(tmpVertexA);
      // fprintf(stderr, " new loop\n");
      for (int vertexIndex = 1; vertexIndex < container->Size(); vertexIndex++) {
        // fprintf(stderr, " . anfang...\n");
        // tmpVertexB = piter.Value();
        tmpVertexB = container->ElementAt(vertexIndex);
        // fprintf(stderr, "%f %f . %f %f\n", tmpVertexA[0], tmpVertexA[1], tmpVertexB[0], tmpVertexB[1]);
        //  const auto indexLineB =
        auto indexLineB = outputImagePtr->TransformPhysicalPointToIndex(tmpVertexB);
        // ContinuousIndexType indexLineB;
        // outputImagePtr->TransformPhysicalPointToContinuousIndex(tmpVertex, indexLineB);
        //  coordinates are indexA - indexB and indexLineA - indexLineB
        _Point2D A, B, C, D;
        if (0) {
          A.x = indexA[0];
          A.y = indexA[1];
          B.x = indexB[0];
          B.y = indexB[1];
          C.x = indexLineA[0];
          C.y = indexLineA[1];
          D.x = indexLineB[0];
          D.y = indexLineB[1];
        } else {
          A.x = floatIndexA[0];
          A.y = floatIndexA[1];
          B.x = floatIndexB[0];
          B.y = floatIndexB[1];
          C.x = tmpVertexA[0];
          C.y = tmpVertexA[1];
          D.x = tmpVertexB[0];
          D.y = tmpVertexB[1];
        }
        /*
        We have a problem here for lines that intersect a polygon vertex. In that case an intersection is
        detected but there are two cases. The line could be detected twice at the same position. In that case
        there are two possibilities. We can either have an ascending line segment or a peak line segment. In
        both cases we will have to be either inside or outside.
        */
        if (doIntersect(A, B, C, D)) {
          if (std::find(lineHitMemory.begin(), lineHitMemory.end(), vertexIndex) != lineHitMemory.end()) {
            // we have hit that polygon side already, only hit once per scanline so ignore that intersection for now
          } else {
            // start flipping the bird (odd flips is inside, even flips is outside)
            outside = !outside;
            //  remember this intersection
            lineHitMemory.push_back(vertexIndex);
          }
        }

        indexLineA[0] = indexLineB[0];
        indexLineA[1] = indexLineB[1];
        tmpVertexA[0] = tmpVertexB[0];
        tmpVertexA[1] = tmpVertexB[1];
      }
      if (outside) {
        imitLine.Set(b_val);
      } else {
        imitLine.Set(f_val);
      }

      indexA[0] = indexB[0];
      indexA[1] = indexB[1];
      floatIndexA[0] = floatIndexB[0];
      floatIndexA[1] = floatIndexB[1];
      ++imitLine;
    }
    imitLine.NextLine();
  }

  // use the line iterator
  /* for (int z = startz; z < endz; z++) {
    for (int y = starty; y < endy; y++) {
      bool outside = true;
      for (int x = startx; x < endx; x++) { // x is the fastest running index
        // two lines, do we cross them?
        piter = container->Begin();
        while (piter != container->End()) {
          tmpVertex = piter.Value();
          const auto tmpIndex = outputImagePtr->TransformPhysicalPointToIndex(tmpVertex);
          // coordinates are
          ++piter;
        }
      }
    }
  }
  */

  /*  // reset piter
    piter = container->Begin();

  bool pflag;

  PixelType zero_val = NumericTraits<PixelType>::ZeroValue();
  auto u_val = static_cast<PixelType>(0);
  auto b_val = static_cast<PixelType>(2);
  auto f_val = static_cast<PixelType>(255);
  outputImagePtr->FillBuffer(u_val);

  pstartVertex = piter.Value();

  tmpVertex = pstartVertex;
  ++piter;

  ImageIndexType tmpImageIndex;
  tmpImageIndex.Fill(0);

  ImageLineIteratorType imit(outputImagePtr, outputImagePtr->GetLargestPossibleRegion());
  imit.SetDirection(0);

  itkDebugMacro(<< "Generating the mask defined by the polyline.....");

  while (piter != container->End()) {
    pflag = false;
    startVertex = tmpVertex;
    endVertex = piter.Value();

    const auto startImageIndex = outputImagePtr->TransformPhysicalPointToIndex(startVertex);
    const auto endImageIndex = outputImagePtr->TransformPhysicalPointToIndex(endVertex);

    // itkDebugMacro(<<"Projection image (index,physical
    // coordinate):"<<startImageIndex<<","<<startVertex<<std::endl);

    if (endImageIndex[1] > startImageIndex[1]) {
      pflag = true;
    }

    LineIteratorType it(outputImagePtr, startImageIndex, endImageIndex);
    it.GoToBegin();

    while (!it.IsAtEnd()) {
      tmpImageIndex[0] = it.GetIndex()[0];
      tmpImageIndex[1] = it.GetIndex()[1];

      // initialize imit using it
      imit.SetIndex(tmpImageIndex);
      while (!imit.IsAtEndOfLine()) {
        if (pflag) {
          if (imit.Get() == u_val) {
            imit.Set(f_val);
          }
        } else {
          imit.Set(b_val);
        }
        ++imit;
      }
      ++it;
    }
    tmpVertex = endVertex;
    ++piter;
  }

  //
  pflag = false;
  startVertex = tmpVertex;
  endVertex = pstartVertex;

  const auto startImageIndex = outputImagePtr->TransformPhysicalPointToIndex(startVertex);
  const auto endImageIndex = outputImagePtr->TransformPhysicalPointToIndex(endVertex);

  if (endImageIndex[1] > startImageIndex[1]) {
    pflag = true;
  }

  LineIteratorType it(outputImagePtr, startImageIndex, endImageIndex);
  it.GoToBegin();

  while (!it.IsAtEnd()) {
    tmpImageIndex[0] = it.GetIndex()[0];
    tmpImageIndex[1] = it.GetIndex()[1];

    // initialize imit using it
    imit.SetIndex(tmpImageIndex);
    while (!imit.IsAtEndOfLine()) {
      if (pflag) {
        if (imit.Get() == u_val) {
          imit.Set(f_val);
        }
      } else {
        imit.Set(b_val);
      }
      ++imit;
    }
    ++it;
  }

  */

  /* Mask the input image with the mask generated */
  //InputImageConstIteratorType inputI(inputImagePtr, inputImagePtr->GetLargestPossibleRegion());
  OutputImageIteratorType outputI(outputImagePtr, outputImagePtr->GetLargestPossibleRegion());
  //inputI.GoToBegin();
  outputI.GoToBegin();
  while (!outputI.IsAtEnd()) {
    if (outputI.Get() == f_val) {
      outputI.Set(1 /*inputI.Get()*/);  // we don't want to mask the pixel values, we want the mask itself (value inside should be 1)
    } else {
      outputI.Set(zero_val);
    }
    //++inputI;
    ++outputI;
  }
}
} // end namespace itk
#endif
