#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#define jbyte int
using namespace std;


void writeImage(std::string fname, vector<jbyte> &img) {
  ofstream file(fname);
  for(int n=0; n<img.size(); n++) file<<(int)img[n]<<",";
  file.close();
}
/*
void writeImage(std::string fname, vector<int> &img) {
  ofstream file(fname);
  for(int n=0; n<img.size(); n++) file<<(int)img[n]<<",";
  file.close();
}
*/

std::vector<double> generateRange(double min_val, double max_val, double step) {
    std::vector<double> range;
    int count = static_cast<int>((max_val - min_val) / step);  // Calculate the number of steps
    for (int i = 0; i <= count; ++i) {
        double value = min_val + i * step;  // Calculate the floating-point value
        if (value <= max_val) {
            range.push_back(value);
        }
    }
    return range;
}
std::vector<double> concatenate(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> result(vec1.size() + vec2.size());
    std::copy(vec1.begin(), vec1.end(), result.begin());
    std::copy(vec2.begin(), vec2.end(), result.begin() + static_cast<int>(vec1.size()));
    return result;
}
std::vector<double> houghTransform(jbyte* cannyData, int width, int height, int rho, double theta,
                                   double bins1[], double bins2[], double thresh) {
    // Bin Format - {min_rho, max_rho, min_theta, max_theta}
    std::vector<double> thetaAngles1 = generateRange(bins1[2], bins1[3], theta);
    std::vector<double> thetaAngles2 = generateRange(bins2[2], bins2[3], theta);
    std::vector<double> thetaAngles = concatenate(thetaAngles1, thetaAngles2);

    std::vector<double> rhoValues1 = generateRange(bins1[0], bins1[1], rho);
    std::vector<double> rhoValues2 = generateRange(bins2[0], bins2[1], rho);
    std::vector<double> rhoValues = concatenate(rhoValues1, rhoValues2);

    std::vector<int> accumulator(thetaAngles.size() * rhoValues.size(), 0);

    std::vector<double> cosValues;
    std::vector<double> sinValues;
    for (auto& t: thetaAngles) {
        cosValues.push_back(std::cos(t));
        sinValues.push_back(std::sin(t));
    }
    
    
    cout << "#rho="<<rhoValues.size()<<", #theta="<<thetaAngles.size()<<endl;
    /*
    for(int n=0; n<rhoValues.size(); n++) cout<<rhoValues[n]<<", ";
    cout<<endl;
    for(int n=0; n<thetaAngles.size(); n++) cout<<thetaAngles[n]<<", ";
    cout<<endl;
    */
    int count = 0, maxCount = 100;
    for (int y = 1; y < height-1; ++y) {
//        if (count > maxCount) break;
        for (int x = 1; x < width-1; ++x) {
            if (0 < cannyData[y * width + x]) {
//                if (count > maxCount) break;
                for (int t = 0; t < thetaAngles.size(); ++t) {
                    double currentRho = y * cosValues[t] + x * sinValues[t];

                    int r = -1;
                    if ((bins2[0] <= currentRho) && (currentRho<= bins2[1])) {
                      r = round((currentRho-bins2[0])/rho) + rhoValues1.size();
                    }
                    if ((bins1[0] <= currentRho) && (currentRho<= bins1[1])) {
                      r = round((currentRho-bins1[0])/rho);
                    }
                    if(r<0) continue;
                        /*
                        int r = 0;
                        double minDiff = std::numeric_limits<double>::max();

                        for (int i = 0; i < rhoValues.size(); ++i) {
                            double rhoDiff = std::abs(currentRho - rhoValues[i]);
                            if (rhoDiff < minDiff) {
                                minDiff = rhoDiff;
                                r = i;
                            }
                        }
                        */
                        accumulator[r * thetaAngles.size() + t] += 1;
                    }
                    
                }
                count++;
            }
        }
    
    cout<<"count="<<count<<endl;
    writeImage("acc.csv", accumulator);

    std::vector<double> finalRho, finalTheta;
    int minVal = *std::min_element(accumulator.begin(), accumulator.end());
    int maxVal = *std::max_element(accumulator.begin(), accumulator.end());
    cout << "Max="<<maxVal<<", Min="<<minVal<<endl;
    for (int r = 0; r < rhoValues.size(); ++r) {
        for (int t = 0; t < thetaAngles.size(); ++t) {
            int acc = accumulator[r * thetaAngles.size() + t];
            int accS = (acc - minVal) / (maxVal - minVal);

            if (thresh <= accS) {
                finalRho.push_back(rhoValues[r]);
                finalTheta.push_back(thetaAngles[t]);
                cout << "R="<<rhoValues[r]<<", T="<<thetaAngles[t]<<", Acc="<<acc<<endl;
            }
        }
    }

//    std::vector<double> concatenated = concatenate(finalRho, finalTheta);
    std::vector<double> concatenated(finalRho.size() + finalTheta.size());
    std::copy(finalRho.begin(), finalRho.end(), concatenated.begin());
    std::copy(finalTheta.begin(), finalTheta.end(), concatenated.begin() + static_cast<int>(finalRho.size()));

    return concatenated;
}

void thresholding(const uint8_t* nmsData, jbyte* threshData, int width, int height,
                  int weakPixel, int strongPixel, double lowThreshRatio, double highThreshRatio) {

    uint8_t max = std::numeric_limits<uint8_t>::max();
    for (int i = 0; i < width * height; ++i)
        if (max < nmsData[i]) max = nmsData[i];

    auto highThresh = static_cast<uint8_t>(max * highThreshRatio);
    auto lowThresh = static_cast<uint8_t>(highThresh * lowThreshRatio);

    for (int y = 1; y < height-1; ++y) {
        for (int x = 1; x < width-1; ++x) {
            if (highThresh <= nmsData[y * width + x])
                threshData[y * width + x] = static_cast<jbyte>(strongPixel);
            else if (lowThresh <= nmsData[y * width + x] &&
                    nmsData[y * width + x] < highThresh)
                threshData[y * width + x] = static_cast<jbyte>(weakPixel);
            else threshData[y * width + x] = static_cast<jbyte>(0);
        }
    }
}
void convolve2D(const jbyte* nv21, jbyte* res, int width, int height, int** kernel, const int scale=1) {
    for (int y = 0; y < height; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int sum = 0;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    sum += kernel[i + 1][j + 1] * nv21[(y + i) * width + (x + j)];
                }
            }
            res[y * width + x] = static_cast<jbyte>(sum/scale);
        }
    }
}
void nonMaxSuppression(const jbyte* magniData, const jbyte* thetaData, jbyte* nmsData, int width, int height) {
    jbyte q, r;
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            try {
                q = static_cast<jbyte>(255);
                r = static_cast<jbyte>(255);

                jbyte angle = thetaData[y * width + x];
                jbyte magni = magniData[y * width + x];

                // Angle 0
                if( (0 <= angle && angle < 22.5) || (157.5 <= angle && angle <= 180.0)) {
                    q = magniData[(y + 1) * width + x];
                    r = magniData[(y - 1) * width + x];
                }

                // Angle 45
                if (22.5 <= angle && angle < 67.5) {
                    q = magniData[(y - 1) * width + (x + 1)];
                    r = magniData[(y + 1) * width + (x - 1)];
                }

                // Angle 90
                if (67.5 <= angle && angle < 112.5) {
                    q = magniData[y * width + (x + 1)];
                    r = magniData[y * width + (x - 1)];
                }

                // Angle 135
                if (112.5 <= angle && angle < 157.5) {
                    q = magniData[(y - 1) * width + (x - 1)];
                    r = magniData[(y + 1) * width + (x + 1)];
                }

                if (q <= magni && r <= magni) nmsData[y * width + x] = magni;
                else nmsData[y * width + x] = static_cast<jbyte>(0);
            }
            catch (const std::exception &e) {
                // Pass
            }
        }
    }
}

void convolve1DRow(const jbyte* nv21, jbyte* res, int width, int height, const int kernel[], const int scale=1) {
    for (int y = 0; y < height; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int sum = 0;
            for (int i = -1; i <= 1; ++i) {
                sum += kernel[i + 1] * nv21[y * width + (x + i)];
            }
            res[y * width + x] = static_cast<jbyte>(sum/scale);
        }
    }
}
void convolve1DCol(const jbyte* nv21, jbyte* res, int width, int height, const int kernel[], const int scale=1) {
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 0; x < width; ++x) {
            int sum = 0;
            for (int i = -1; i <= 1; ++i) {
                sum += kernel[i + 1] * nv21[(y + i) * width + x];
            }
            res[y * width + x] = static_cast<jbyte>(sum/scale);
        }
    }
}
void hysteresis(const uint8_t* threshData, jbyte* cannyData, int width, int height,  int weakPixel, int strongPixel) {
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            if (threshData[y * width + x] == weakPixel) {
                try {
                    if ((threshData[(y - 1) * width + x] == strongPixel) ||
                        (threshData[(y + 1) * width + x] == strongPixel) ||

                        (threshData[y * width + (x + 1)] == strongPixel) ||
                        (threshData[(y - 1) * width + (x + 1)] == strongPixel) ||
                        (threshData[(y + 1) * width + (x + 1)] == strongPixel) ||

                        (threshData[y * width + (x - 1)] == strongPixel) ||
                        (threshData[(y - 1) * width + (x - 1)] == strongPixel) ||
                        (threshData[(y + 1) * width + (x - 1)] == strongPixel))
                        cannyData[y * width + x] = static_cast<jbyte>(strongPixel);
                    else cannyData[y * width + x] = static_cast<jbyte>(0);
                }
                catch (const std::exception &e) {
                    // Pass
                }
            }
            else cannyData[y * width + x] = static_cast<jbyte>(threshData[y * width + x]);
        }
    }
}

std::vector<std::string> adv_tokenizer(string s, char del) {
  std::vector<string> words;
  stringstream ss(s);
  string word;
  while (!ss.eof()) {
      getline(ss, word, del);
      words.push_back(word);
  }
  return words;
}

void getNumbers(std::string line, std::vector<jbyte> &tmp) {
  std::vector<string> words = adv_tokenizer(line, ',');
  for(int n=0; n<words.size(); n++) {
    int val;
    sscanf(words[n].c_str(), "%d", &val);
    unsigned char c = (unsigned char)val;
    // cout<<words[n].c_str()<<", "<<val<<", "<<(int)c<<endl;
    tmp.push_back((jbyte)val);
  }
}

vector<jbyte> readImage(string fname) {
  ifstream file(fname);
  int nlines = 0;
  std::vector<jbyte> tmp;
  if(file.is_open()) {
    string line;
    while(getline(file, line)) {
      getNumbers(line, tmp);
    }
  }
  return tmp;
}


void writeImage(std::string fname, jbyte* img, int W, int H) {
  ofstream file(fname);
  for(int n=0; n<W*H; n++) file<<(int)img[n]<<",";
  file.close();
}


int main(int argc, char* argv[]) {
  /*
  std::string s("1,2,4,255");
  cout<<s.c_str()<<endl;
  std::vector<jbyte> tmp;
  getNumbers(s, tmp);
  for(int n=0; n<tmp.size(); n++) {
    cout<<(int)tmp[n]<<";";
  }
  cout<<endl;
  */
  
  vector<jbyte> tmp = readImage(string("exTest_80x80.csv"));
  
  jbyte* nv21 = tmp.data();
  
  for(int n=0; n<10; n++) cout<<(int)nv21[n]<<",";
  cout<<endl;
  for(int n=0; n<10; n++) cout<<(int)nv21[n+100]<<",";
  cout<<endl;
  for(int n=0; n<10; n++) cout<<(int)nv21[n+1000]<<",";
  cout<<endl;
  
  int derivative[3] = {-1, 0, 1};        
  int gauss_filter[3] = {1, 2, 1};
  int width=80, height=80;
  int N = width*height;
  jbyte *nv21TempData = new jbyte[N];
  jbyte *nv21GxData = new jbyte[N];
  jbyte *nv21GyData = new jbyte[N];
  jbyte *nv21MagniData = new jbyte[N];
  jbyte *nv21ThetaData = new jbyte[N];
  jbyte *nv21NMSData = new jbyte[N];
  for(int n=0; n<width*height; n++) {
    nv21TempData[n] = 0;
    nv21GxData[n] = 0;
    nv21GyData[n] = 0;
    nv21MagniData[n] = 0;
    nv21ThetaData[n] = 0;
    nv21NMSData[n] = 0;
  }
  
  convolve1DRow(nv21, nv21TempData,  width, height, derivative, 1);
  convolve1DCol(nv21TempData, nv21GxData, width, height, gauss_filter, 4);

  convolve1DRow(nv21, nv21TempData,  width, height, gauss_filter, 4);
  convolve1DCol(nv21TempData, nv21GyData, width, height, derivative, 1);
  
  for (int y = 1; y < height - 1; ++y) {
    for (int x = 1; x < width - 1; ++x) {
        double mag = sqrt(pow(nv21GxData[y * width + x], 2) + pow(nv21GyData[y * width + x], 2));
    
        nv21MagniData[y * width + x] = static_cast<jbyte>((int)mag);
                                                         
        double theta = atan2(nv21GxData[y * width + x], nv21GyData[y * width + x]) * 180.0/M_PI;
        nv21ThetaData[y * width + x] = static_cast<jbyte>((theta < 0) ? theta + 180 : theta);
    }
  }
  
  nonMaxSuppression(nv21MagniData, nv21ThetaData, nv21NMSData, width, height);

  // thresholding(nms, nv21ThreshData, width, height, 40, 100, low_thresh_ratio, high_thresh_ratio);
  
  jbyte* edgeImage = nv21NMSData;
  for(int n=0; n<width*height; n++) if(edgeImage[n]>20) edgeImage[n]=255; else edgeImage[n]=0;
  // for(int n=0; n<width*height; n++) nv21MagniData[n] = 0;
  // for(int n=0; n<height; n++) nv21MagniData[n*width+40] = 200;
  
  double bins1[] = {-46.0, 28.0, -1.67, -0.23};
  double bins2[] = {28.0, 46.0, 0.23, 1.67};
  std::vector<double> HT = houghTransform(edgeImage, width, height, 2, 5*3.14/180.0, bins1, bins2, 0.60);
  
  cout << endl;
  for(int n=0; n<HT.size(); n++) {
    cout<<HT[n]<<", ";
  }
  cout << endl;
  
  writeImage("tmp.csv", edgeImage, width, height);
  
  return 0;
}

/*
extern "C"
JNIEXPORT jobjectArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_toSobelConv1D(
        JNIEnv *env, jclass clazz, jbyteArray nv21Array, jint width, jint height) {
    try {
        jbyte *nv21 = env->GetByteArrayElements(nv21Array, nullptr);

        jbyteArray nv21Gx = env->NewByteArray(width * height);
        jbyteArray nv21Gy = env->NewByteArray(width * height);
        jbyteArray nv21Temp = env->NewByteArray(width * height);
        jbyteArray nv21Magni = env->NewByteArray(width * height * 3/2);
        jbyteArray nv21Theta = env->NewByteArray(width * height * 3/2);

        jbyte *nv21GxData = env->GetByteArrayElements(nv21Gx, nullptr);
        jbyte *nv21GyData = env->GetByteArrayElements(nv21Gy, nullptr);
        jbyte *nv21TempData = env->GetByteArrayElements(nv21Temp, nullptr);
        jbyte *nv21MagniData = env->GetByteArrayElements(nv21Magni, nullptr);
        jbyte *nv21ThetaData = env->GetByteArrayElements(nv21Theta, nullptr);

        jobjectArray resultArray = env->NewObjectArray(2,
                                                       env->FindClass("[B"), nullptr);

        std::memset(nv21MagniData, 128, width * height * 3/2);
        std::memset(nv21ThetaData, 128, width * height * 3/2);

        // Separable 1D Sobel Kernels
        int derivative[3] = {-1, 0, 1};
        int gauss_filter[3] = {1, 2, 1};

        convolve1DRow(nv21, nv21TempData,  width, height, derivative);
        convolve1DCol(nv21TempData, nv21GxData, width, height, gauss_filter);

        convolve1DRow(nv21, nv21TempData,  width, height, gauss_filter);
        convolve1DCol(nv21TempData, nv21GyData, width, height, derivative);

        // Calculate magnitude (G) and theta (direction)
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                nv21MagniData[y * width + x] = static_cast<jbyte>(sqrt(pow(nv21GxData[y * width + x], 2) +
                                                                       pow(nv21GyData[y * width + x], 2)));
                double theta = atan2(nv21GxData[y * width + x], nv21GyData[y * width + x]) * 180.0/M_PI;
                nv21ThetaData[y * width + x] = static_cast<jbyte>((theta < 0) ? theta + 180 : theta);
            }
        }

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21Magni, nv21MagniData, 0);
        env->ReleaseByteArrayElements(nv21Theta, nv21ThetaData, 0);

        env->SetObjectArrayElement(resultArray, 0, nv21Magni);
        env->SetObjectArrayElement(resultArray, 1, nv21Theta);

        return resultArray;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jobjectArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_toSobelConv2D(
        JNIEnv *env, jclass clazz, jbyteArray nv21Array, jint width, jint height) {
    try {
        jbyte *nv21 = env->GetByteArrayElements(nv21Array, nullptr);

        jbyteArray nv21Magni = env->NewByteArray(width * height * 3/2);
        jbyteArray nv21Theta = env->NewByteArray(width * height * 3/2);

        jbyte *nv21MagniData = env->GetByteArrayElements(nv21Magni, nullptr);
        jbyte *nv21ThetaData = env->GetByteArrayElements(nv21Theta, nullptr);

        jobjectArray resultArray = env->NewObjectArray(2,
                                                       env->FindClass("[B"), nullptr);

        std::memset(nv21MagniData, 128, width * height * 3/2);
        std::memset(nv21ThetaData, 128, width * height * 3/2);

        // Sobel kx-Ky kernels
        int GxKernel[3][3] = {
                {-1, 0, 1},
                {-2, 0, 2},
                {-1, 0, 1}
        };

        int GyKernel[3][3] = {
                {-1, -2, -1},
                { 0,  0,  0},
                { 1,  2,  1}
        };

        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                int sumX = 0, sumY = 0;

                // Convolution operation with Gx and Gy kernels
                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        sumX += GxKernel[i + 1][j + 1] * nv21[(y + i) * width + (x + j)];
                        sumY += GyKernel[i + 1][j + 1] * nv21[(y + i) * width + (x + j)];
                    }
                }

                // Calculation for gradient magnitude (G) and gradient direction (theta)
                double theta = atan2(sumY, sumX) * 180.0/M_PI;
                nv21ThetaData[y * width + x] = static_cast<jbyte>((theta < 0) ? theta + 180 : theta);
                nv21MagniData[y * width + x] = static_cast<jbyte>(sqrt(sumX * sumX + sumY * sumY));
            }
        }

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21Magni, nv21MagniData, 0);
        env->ReleaseByteArrayElements(nv21Theta, nv21ThetaData, 0);

        env->SetObjectArrayElement(resultArray, 0, nv21Magni);
        env->SetObjectArrayElement(resultArray, 1, nv21Theta);

        return resultArray;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_toGaussianFilterConv1D(
        JNIEnv *env, jclass clazz, jbyteArray nv21Array, jint width, jint height) {
    try {
        jbyte *nv21 = env->GetByteArrayElements(nv21Array, nullptr);

        jbyteArray nv21Temp = env->NewByteArray(width * height);
        jbyteArray nv21Grad = env->NewByteArray(width * height * 3 / 2);

        jbyte *nv21TempData = env->GetByteArrayElements(nv21Temp, nullptr);
        jbyte *nv21GradData = env->GetByteArrayElements(nv21Grad, nullptr);

        std::memset(nv21GradData, 128, width * height * 3/2);

        int gauss_filter[5] = {1, 4, 6, 4, 1};
        convolve1DCol(nv21, nv21TempData, width, height, gauss_filter, 16);
        convolve1DRow(nv21TempData, nv21GradData,  width, height, gauss_filter, 16);

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21Grad, nv21GradData, 0);

        return nv21Grad;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_toGaussianFilterConv2D(
        JNIEnv *env, jclass clazz, jbyteArray nv21Array, jint width, jint height) {
    try {
        jbyte *nv21 = env->GetByteArrayElements(nv21Array, nullptr);

        jbyteArray nv21Gauss = env->NewByteArray(width * height * 3/2);

        jbyte *nv21GaussData = env->GetByteArrayElements(nv21Gauss, nullptr);
        std::memset(nv21GaussData + (width * height), 128, width * height * 3/2);

        // Allocate memory for a 3x3 matrix
        int rows = 5;
        int cols = 5;
        auto** kernel = new int*[rows];
        for (int i = 0; i < rows; ++i)
            kernel[i] = new int[cols];

        // Initialize the kernel with the given values
        int values[5][5] = {{1, 4,  6,  4,  1},
                            {4, 16, 24, 16, 4},
                            {6, 24, 36, 24, 6},
                            {4, 16, 24, 16, 4},
                            {1, 4,  6,  4,  1}};

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                kernel[i][j] = values[i][j];
            }
        }

        convolve2D(nv21, nv21GaussData, width, height, kernel, 256);

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21Gauss, nv21GaussData, 0);

        return nv21Gauss;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_nonMaxSuppression(
        JNIEnv *env, jclass clazz, jobjectArray sobelArray, jint width, jint height) {
    try {
        auto magniArray = (jbyteArray) env->GetObjectArrayElement(sobelArray, 0);
        auto thetaArray = (jbyteArray) env->GetObjectArrayElement(sobelArray, 0);

        jbyte *nv21MagniData = env->GetByteArrayElements(magniArray, nullptr);
        jbyte *nv21ThetaData = env->GetByteArrayElements(thetaArray, nullptr);

        jbyteArray nv21NMS = env->NewByteArray(width * height * 3/2);

        jbyte *nv21NMSData = env->GetByteArrayElements(nv21NMS, nullptr);

        std::memset(nv21NMSData, 128, width * height * 3/2);

        nonMaxSuppression(nv21MagniData, nv21ThetaData, nv21NMSData, width, height);

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21NMS, nv21NMSData, 0);

        return nv21NMS;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_thresholding(
        JNIEnv *env, jclass clazz, jbyteArray nmsArray, jint width, jint height, jint weak_pixel,
        jint strong_pixel, jdouble low_thresh_ratio, jdouble high_thresh_ratio) {
    try {
//        jbyte *nms = env->GetByteArrayElements(nmsArray, nullptr);
        int len = width * height;
        auto* nms = new uint8_t[len];
        env->GetByteArrayRegion (nmsArray, 0, len, reinterpret_cast<jbyte*>(nms));

        jbyteArray nv21Thresh = env->NewByteArray(width * height * 3/2);

        jbyte *nv21ThreshData = env->GetByteArrayElements(nv21Thresh, nullptr);

        std::memset(nv21ThreshData, 128, width * height * 3/2);

        thresholding(nms, nv21ThreshData, width, height, weak_pixel,
                     strong_pixel, low_thresh_ratio, high_thresh_ratio);

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21Thresh, nv21ThreshData, 0);

        return nv21Thresh;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_hysteresis(
        JNIEnv *env, jclass clazz, jbyteArray threshArray, jint width, jint height, jint weak_pixel,
        jint strong_pixel) {
    try {
//        jbyte *thresh = env->GetByteArrayElements(threshArray, nullptr);
        int len = width * height;
        auto* thresh = new uint8_t[len];
        env->GetByteArrayRegion (threshArray, 0, len, reinterpret_cast<jbyte*>(thresh));

        jbyteArray nv21Canny = env->NewByteArray(width * height * 3/2);

        jbyte *nv21CannyData = env->GetByteArrayElements(nv21Canny, nullptr);

        std::memset(nv21CannyData, 128, width * height * 3/2);

        hysteresis(thresh, nv21CannyData, width, height, weak_pixel, strong_pixel);

        // Release the nv21 array
        env->ReleaseByteArrayElements(nv21Canny, nv21CannyData, 0);

        return nv21Canny;
//

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}

extern "C"
JNIEXPORT jdoubleArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_houghTransform(
        JNIEnv *env, jclass clazz, jbyteArray cannyArray, jint width, jint height, jint rho, jdouble theta,
        jdoubleArray bins1Array, jdoubleArray bins2Array, jdouble thresh) {
    try {
//        jbyte *canny = env->GetByteArrayElements(cannyArray, nullptr);
        int len = width * height;
        auto* canny = new uint8_t[len];
        env->GetByteArrayRegion (cannyArray, 0, len, reinterpret_cast<jbyte*>(canny));

        jdouble *bins1 = env->GetDoubleArrayElements(bins1Array, nullptr);
        jdouble *bins2 = env->GetDoubleArrayElements(bins2Array, nullptr);

        std::vector<double> houghTransformData = houghTransform(canny, width, height, rho, theta,
                                                                bins1, bins2, thresh);

        int transformedArraySize = static_cast<int>(houghTransformData.size());
        jdoubleArray transformedArray = env->NewDoubleArray(transformedArraySize);
        env->SetDoubleArrayRegion(transformedArray, 0,transformedArraySize, houghTransformData.data());
        return transformedArray;

    } catch (const std::exception &e) {
        jclass exceptionClass = env->FindClass("java/lang/Exception");
        env->ThrowNew(exceptionClass, e.what());
    }

    return nullptr;
}
*/