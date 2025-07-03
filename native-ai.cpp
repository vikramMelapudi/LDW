// #include <jni.h>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <fstream>
// #include <android/log.h>

#define intArray std::vector<int>
#define floatArray std::vector<double>
// #define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, "LDW", __VA_ARGS__)
// #define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, "LDW", __VA_ARGS__)

using namespace std;

void nonMaxSuppression(intArray &magni, intArray &theta, intArray &nmsData, int width, int height) {
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            try {
                int q = 255, r = 255;

                int mag = magni[y * width + x];
                int angle = theta[y * width + x];
                
                // Angle 0
                if( (0 <= angle && angle < 22.5) || (157.5 <= angle && angle <= 180.0)) {
                    q = magni[y* width + x+1];
                    r = magni[y* width + x-1];
                }

                // Angle 45
                if (22.5 <= angle && angle < 67.5) {
                    q = magni[(y + 1) * width + (x - 1)];
                    r = magni[(y - 1) * width + (x + 1)];
                }

                // Angle 90
                if (67.5 <= angle && angle < 112.5) {
                    q = magni[(y+1) * width + x];
                    r = magni[(y-1) * width + x];
                }

                // Angle 135
                if (112.5 <= angle && angle < 157.5) {
                    q = magni[(y - 1) * width + (x - 1)];
                    r = magni[(y + 1) * width + (x + 1)];
                }

                if (q <= mag && r <= mag) nmsData[y * width + x] = mag;
                else nmsData[y * width + x] = 0;
            }
            catch (const std::exception &e) {
                // Pass
            }
        }
    }
}
void convolve2D(intArray &nv21, intArray &res, int width, int height, int** kernel, int size=3, int scale=1) {
    int half = size/2;
    for (int y = half; y < height - half; ++y) {
        for (int x = half; x < width - half; ++x) {
            int sum = 0;
            for (int ky = -half; ky <= half; ++ky) {
                for (int kx = -half; kx <= half; ++kx) {
                    int pixel = nv21[(y + ky) * width + (x + kx)];
                    sum += pixel * kernel[ky + half][kx + half];
                }
            }
            res[y * width + x] = sum/scale;
        }
    }
}
void convolve1DRow(intArray &nv21, intArray &res, int width, int height, int* kernel, int size=3, int scale=1) {
    int half = size/2;
    for (int y = 0; y < height; ++y) {
        for (int x = half; x < width - half; ++x) {
            int sum = 0;
            for (int kx = -half; kx <= half; ++kx) {
                int pixel = nv21[y * width + (x + kx)];
                sum += pixel * kernel[kx + half];
            }
            res[y * width + x] = sum/scale;
        }
    }
}
void convolve1DCol(intArray &nv21, intArray &res, int width, int height, int* kernel, int size=3, int scale=1) {
    int half = size/2;
    for (int y = half; y < height - half; ++y) {
        for (int x = 0; x < width; ++x) {
            int sum = 0;
            for (int ky = -half; ky <= half; ++ky) {
                int pixel = nv21[(y + ky) * width + x];
                sum += pixel * kernel[ky + half];
            }
            res[y * width + x] = sum/scale;
        }
    }
}
void hysteresis(intArray &threshData, intArray &cannyData, int width, int height,  int weakPixel, int strongPixel) {
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
                        cannyData[y * width + x] = strongPixel;
                    else cannyData[y * width + x] = 0;
                }
                catch (const std::exception &e) {
                    // Pass
                }
            }
            else cannyData[y * width + x] = threshData[y * width + x];
        }
    }
}
void convolve2DSobel(intArray &nv21, intArray &magni, intArray &theta, int width, int height, int** gX, int** gY, int size=3, int scale=1) {
    int half = size/2;
    for (int y = half; y < height - half; ++y) {
        for (int x = half; x < width - half; ++x) {
            int sumX = 0, sumY = 0;
            for (int ky = -half; ky <= half; ++ky) {
                for (int kx = -half; kx <= half; ++kx) {
                    int pixel = nv21[(y + ky) * width + (x + kx)];
                    sumX += pixel * gX[ky + half][kx + half]/scale;
                    sumX += pixel * gY[ky + half][kx + half]/scale;
                }
            }
            double t = atan2(sumY, sumX) * 180.0/M_PI;
            theta[y * width + x] = static_cast<int>((t < 0) ? t + 180 : t);
            magni[y * width + x] = static_cast<int>(sqrt(pow(sumX, 2) +
                                                            pow(sumY, 2)));
        }
    }
}
void convolve1DSobel(intArray &nv21, intArray &magni, intArray &theta, int width, int height, int* derivative, int* gauss_filter, int size=3, int scale=1) {
    int length = width * height;
    intArray temp(length, 0);
    intArray gX(length), gY(length);

    temp = intArray(length, 0);
    convolve1DCol(nv21, temp,   width, height, gauss_filter, size, 4);
    convolve1DRow(temp, gX, width, height, derivative,  size, 1);

    // temp = intArray(length, 0);
    if(0) {
      convolve1DRow(nv21, temp,  width, height, gauss_filter, size, 4);
      convolve1DCol(temp, gY,    width, height, derivative,  size, 1);
    } else {
      convolve1DCol(nv21, temp,    width, height, derivative,  size, 1);
      convolve1DRow(temp, gY,  width, height, gauss_filter, size, 4);
    }
    temp.clear();

    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int sumX = gX[y * width + x];
            int sumY = gY[y * width + x];
            double t = atan2(sumY, sumX) * 180.0/M_PI;
            int ang = (int) t;
            if(ang<0) ang += 180;
            // theta[y * width + x] = static_cast<int>((t < 0) ? t + 180 : t);
            theta[y * width + x] = ang;
            magni[y * width + x] = static_cast<int>(sqrt(pow(sumX, 2) +
                                                            pow(sumY, 2)));
        }
    }
    gX.clear(), gY.clear();
}
void thresholding(intArray &nmsData, intArray &threshData, int width, int height, int weakPixel, int strongPixel, double lowThreshRatio, double highThreshRatio) {
    int max = *std::max_element(nmsData.begin(), nmsData.end());
    int highThresh = static_cast<int>(max * highThreshRatio + 1);
    int lowThresh = static_cast<int>(highThresh * lowThreshRatio + 1);

    for (int y = 1; y < height-1; ++y) {
        for (int x = 1; x < width-1; ++x) {
            if (highThresh < nmsData[y * width + x])
                threshData[y * width + x] = strongPixel;
            else if (lowThresh < nmsData[y * width + x] &&
                     nmsData[y * width + x] < highThresh)
                threshData[y * width + x] = weakPixel;
            else threshData[y * width + x] = 0;
        }
    }
}

floatArray concatenate(const floatArray& vec1, const floatArray& vec2) {
    floatArray result(vec1.size() + vec2.size());
    std::copy(vec1.begin(), vec1.end(), result.begin());
    std::copy(vec2.begin(), vec2.end(), result.begin() + static_cast<int>(vec1.size()));
    return result;
}
floatArray generateRange(double min_val, double max_val, double step) {
    floatArray range;
    int count = static_cast<int>((max_val - min_val) / step);  // Calculate the number of steps
    for (int i = 0; i <= count; ++i) {
        double value = min_val + i * step;  // Calculate the floating-point value
        if (value <= max_val) {
            range.push_back(value);
        }
    }
    return range;
}
floatArray houghTransform(intArray &cannyData, int width, int height, int rho, double theta,  double bins1[], double bins2[], double thresh) {
    // Bin Format - {min_rho, max_rho, min_theta, max_theta}
    floatArray thetaAngles1 = generateRange(bins1[2], bins1[3], theta);
    floatArray thetaAngles2 = generateRange(bins2[2], bins2[3], theta);
    floatArray thetaAngles = concatenate(thetaAngles1, thetaAngles2);

    floatArray rhoValues1 = generateRange(bins1[0], bins1[1], rho);
    floatArray rhoValues2 = generateRange(bins2[0], bins2[1], rho);
    floatArray rhoValues = concatenate(rhoValues1, rhoValues2);
    
    /* // Added for debugging bins, and (rho,theta) counts
    for (auto& t: thetaAngles) cout<<t<<",";
    cout<<endl;
    for (auto& t: rhoValues) cout<<t<<",";
    cout<<endl;
    cout<<"#rho="<<rhoValues.size()<<", #Theta="<<thetaAngles.size()<<endl;
    */

    std::vector<int> accumulator(thetaAngles.size() * rhoValues.size(), 0);

    floatArray cosValues;
    floatArray sinValues;
    for (auto& t: thetaAngles) {
        cosValues.push_back(std::cos(t));
        sinValues.push_back(std::sin(t));
    }

//    int count = 0, maxCount = 100;
    for (int y = 1; y < height-1; ++y) {
//        if (count > maxCount) break;
        if (y<20) continue; // omit top rows : sim trapz
        for (int x = 1; x < width-1; ++x) {
            if((x<40)&&(y*0.6<=(40-x))) continue; // omit top left triangle : sim trapz
            if((x>40)&&(y*0.6<=(x-40))) continue; // omit top right triangle : sim trapz
            
            if (0 < cannyData[y * width + x]) {
//                if (count > maxCount) break;
                for (int t = 0; t < thetaAngles.size(); ++t) {
                    double currentRho = y * cosValues[t] + x * sinValues[t];
                    
                    int r = -1;
                    // enforce only correct theta, rho bin is filled
                    if ((thetaAngles[t]>0)&&(bins2[0] <= currentRho) && (currentRho<= bins2[1])) {
                      r = round((currentRho-bins2[0])/rho) + rhoValues1.size();
                    }
                    // enforce only correct theta, rho bin is filled
                    if ((thetaAngles[t]<0)&&(bins1[0] <= currentRho) && (currentRho<= bins1[1])) {
                      r = round((currentRho-bins1[0])/rho);
                    }
                    
                    if(r<0) continue; // no bin to be filled
                    accumulator[r * thetaAngles.size() + t] += 1;
                    
                    /*
                    if ((bins2[0] <= currentRho && currentRho <= bins2[1]) ||
                        (bins1[0] <= currentRho && currentRho <= bins1[1])) {
                        int r = 0;
                        double minDiff = std::numeric_limits<double>::max();

                        for (int i = 0; i < rhoValues.size(); ++i) {
                            double rhoDiff = std::abs(currentRho - rhoValues[i]);
                            if (rhoDiff < minDiff) {
                                minDiff = rhoDiff;
                                r = i;
                            }
                        }

                        accumulator[r * thetaAngles.size() + t] += 1;
                    }
                    */
                    
                }
//                count++;
            }
        }
    }

    floatArray finalAcc, finalRho, finalTheta;
    int minVal = *std::min_element(accumulator.begin(), accumulator.end());
    int maxVal = *std::max_element(accumulator.begin(), accumulator.end());
//    LOGD("Accumulator: Min - %d, Max - %d", minVal, maxVal);

    for (int r = 0; r < rhoValues.size(); ++r) {
        for (int t = 0; t < thetaAngles.size(); ++t) {
            
            int acc = accumulator[r * thetaAngles.size() + t];
            float accScaled = ((float)acc - minVal) / (maxVal - minVal);
            
            if (thresh <= accScaled) {
                finalAcc.push_back(acc);
                finalRho.push_back(rhoValues[r]);
                finalTheta.push_back(thetaAngles[t]);
//                LOGD("[ Row - %.2lf, Theta - %.2lf, Acc - %d ]", rhoValues[r], thetaAngles[t], acc);
            }
        }
    }
    
    // sort (by acc, descending) and get the top 10 entries
    vector<int> index;
    for(int n=0; n<finalAcc.size(); n++) index.push_back(n);
    stable_sort(index.begin(), index.end(), [&finalAcc](size_t i1, size_t i2) {return finalAcc[i1] > finalAcc[i2];});
    float rad2Ang = 180.0/3.147;
    floatArray nDetValues(30, 0);
    for(int n=0; n<index.size(); n++) {
      if(n>=10) break;
      int sn = index[n];
      nDetValues[n*3] = finalAcc[sn];
      nDetValues[n*3+1] = finalRho[sn];
      nDetValues[n*3+2] = finalTheta[sn]*rad2Ang;
    }
    
    std::stringstream ss;
    ss << "[" << minVal << "," << maxVal << "," << finalAcc.size() << ",";
    for (auto& i: nDetValues) ss<<(int)i<<","; //ss << std::fixed << std::setprecision(2) << i << ", ";

    std::string res = ss.str();
    res = res.substr(0, res.size() - 2);
    res += "]";
    
    cout << endl << res << endl;
//      LOGD("%s", res.c_str());

    minVal = *std::min_element(cannyData.begin(), cannyData.end());
    maxVal = *std::max_element(cannyData.begin(), cannyData.end());
//    LOGD("EdgeMap: Min - %d, Max - %d", minVal, maxVal);

    return concatenate(finalRho, finalTheta);
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


void getNumbers(std::string line, std::vector<int> &tmp) {
  std::vector<string> words = adv_tokenizer(line, ',');
  for(int n=0; n<words.size(); n++) {
    int val;
    sscanf(words[n].c_str(), "%d", &val);
    unsigned char c = (unsigned char)val;
    // cout<<words[n].c_str()<<", "<<val<<", "<<(int)c<<endl;
    tmp.push_back((int)val);
  }
}

vector<int> readImage(string fname) {
  ifstream file(fname);
  int nlines = 0;
  std::vector<int> tmp;
  if(file.is_open()) {
    string line;
    while(getline(file, line)) {
      getNumbers(line, tmp);
    }
  }
  return tmp;
}

void writeImage(std::string fname, intArray img, int W, int H) {
  ofstream file(fname);
  for(int n=0; n<W*H; n++) file<<(int)img[n]<<",";
  file.close();
}


int main(int argc, char* argv[]) {
  
  intArray image = readImage(argv[1]);
  cout<<"image size="<<image.size()<<endl;
  int width=80, height=80;
  int length = width*height;
  
  // sobel edge detection
  intArray magniInt(length, 0), thetaInt(length, 0);
  const int kernelSize = 3; 
  if(1) {
    const int valuesG[kernelSize] = {1, 2, 1};
    const int valuesD[kernelSize] = {-1, 0, 1};
    // Perform 1D Convolution
    convolve1DSobel(image, magniInt, thetaInt, width, height, (int *)valuesD, (int *)valuesG, kernelSize, 4);
  
  } else {
    
    const int valuesX[kernelSize][kernelSize] = {
            {-1, 0, 1},
            {-2, 0, 2},
            {-1, 0, 1}
    };
    const int valuesY[kernelSize][kernelSize] = {
            {-1, -2, -1},
            {0,  0,  0},
            {1,  2,  1}
    };
    int **gX = new int *[kernelSize];
    int **gY = new int *[kernelSize];
    for (int i = 0; i < kernelSize; ++i) {
        gX[i] = new int[kernelSize];
        gY[i] = new int[kernelSize];
        for (int j = 0; j < kernelSize; ++j) {
            gX[i][j] = valuesX[i][j];
            gY[i][j] = valuesY[i][j];
        }
    }
    // Perform 2D Convolution
    // convolve2DSobel(image, magniInt, thetaInt, width, height, gX, gY, kernelSize, 4);
    convolve2DSobel(image, magniInt, thetaInt, width, height, (int**)valuesX, (int **)valuesY, kernelSize, 4);
  }
  
  writeImage("out_sobel.dat", magniInt, width, height);
  writeImage("out_sobelAngle.dat", thetaInt, width, height);
  
  /*
  for(int n=0; n<width-1; n++) {
    for(int m=0; m<height-1; m++) {
      magniInt[n+m*width] = max(magniInt[n+m*width],magniInt[n+1+m*width]);
    }
  }
  for(int n=0; n<width-1; n++) {
    for(int m=0; m<height-1; m++) {
      magniInt[n+m*width] = max(magniInt[n+m*width],magniInt[n+(m+1)*width]);
    }
  }
  */
  
  // Perform Non-Max Suppression
  intArray nmsInt(length, 0);
  nonMaxSuppression(magniInt, thetaInt, nmsInt, width, height);
  writeImage("out_nonMax.dat", nmsInt, width, height);
  
  // thresholding
  intArray threshInt(length, 0);
  int weak_pixel=100, strong_pixel=255;
  double low_thresh_ratio=0.09f, high_thresh_ratio=0.17f;
  if (1) {
    thresholding(nmsInt, threshInt, width, height, weak_pixel, strong_pixel,
                     low_thresh_ratio, high_thresh_ratio);
  } else {
    thresholding(magniInt, threshInt, width, height, weak_pixel, strong_pixel,
                     low_thresh_ratio, high_thresh_ratio);
  }
  writeImage("out_thresh.dat", threshInt, width, height);
  
  // hysteresis
  intArray cannyInt(length, 0);
  hysteresis(threshInt, cannyInt, width, height, weak_pixel, strong_pixel);
  writeImage("out_hysteresis.dat", cannyInt, width, height);
  
  intArray &edge = cannyInt;
  
  cout<<"#entries="<<edge.size()<<endl;
  int minVal = *std::min_element(edge.begin(), edge.end());
  int maxVal = *std::max_element(edge.begin(), edge.end());
  cout<<"min="<<minVal<<", max="<<maxVal<<endl;
  
  for(int n=0; n<edge.size(); n++) if(edge[n]>0) edge[n]=255; else edge[n]=0;

  
  int W=80, H=80;
  int rho=4;
  double theta = 5 * 3.14/180;
  // double[] bins1 = {-46.0, -28.0, -1.67, -0.23}, bins2 = {28.0, 46.0, 0.23, 1.67};
  /*
  double bins1[] = {-46.0, -28.0, -1.67, -0.23};
  double bins2[] = {28.0, 46.0, 0.23, 1.67};
  */
  double bins1[] = {-60.0, -5.0, -1.6, -0.94};
  double bins2[] = {20.0, 70.0, 0.94, 1.6};
  
  
  floatArray HT = houghTransform(edge, W, H, rho, theta, bins1, bins2, 0.75);
  
  
  return 0;
}

/*
extern "C"
JNIEXPORT jobjectArray JNICALL
Java_intangles_videotelematics_invision_1face_1detectionandembedding_NativeUtils_toSobelConv1D(
        JNIEnv *env, jclass clazz, jbyteArray nv21Array, jint width, jint height) {
    try {
        jsize length = width * height;
        jbyte *nv21Byte = env->GetByteArrayElements(nv21Array, nullptr);

        // Convert Data to Int32
        intArray nv21Int(length, 0), gXInt(length, 0), gYInt(length, 0);
        for (int i = 0; i < length; ++i) nv21Int[i] = static_cast<int>(nv21Byte[i] & 0xFF);

        // Define Kernel Filter
        const int kernelSize = 3;
        int *derivative = new int[kernelSize];
        int *gauss_filter = new int[kernelSize];
        const int valuesG[kernelSize] = {1, 2, 1};
        const int valuesD[kernelSize] = {-1, 0, 1};
        for (int i= 0; i < kernelSize; ++i) {
            derivative[i] = valuesD[i], gauss_filter[i] = valuesG[i];
        }

        // Perform 1D Convolution
        intArray magniInt(length, 0), thetaInt(length, 0);
        convolve1DSobel(nv21Int, magniInt, thetaInt, width, height, derivative, gauss_filter, kernelSize, 4);

        // Convert Back to UInt8
        jbyteArray magniArray = env->NewByteArray(length * 3 / 2);
        jbyteArray thetaArray = env->NewByteArray(length * 3 / 2);

        jbyte *magniByte = env->GetByteArrayElements(magniArray, nullptr);
        jbyte *thetaByte = env->GetByteArrayElements(thetaArray, nullptr);
        std::memset(magniByte, 128, length * 3 / 2);
        std::memset(thetaByte, 128, length * 3 / 2);

        for (int i = 0; i < length; ++i) {
            magniByte[i] = static_cast<jbyte>(magniInt[i] & 0xFF);
            thetaByte[i] = static_cast<jbyte>(thetaInt[i] & 0xFF);
        }

        // Release Resources
        nv21Int.clear(), magniInt.clear(), thetaInt.clear();
        env->ReleaseByteArrayElements(nv21Array, nv21Byte, 0);
        env->ReleaseByteArrayElements(magniArray, magniByte, 0);
        env->ReleaseByteArrayElements(thetaArray, thetaByte, 0);

        jobjectArray resultArray = env->NewObjectArray(2,
                                                       env->FindClass("[B"), nullptr);
        env->SetObjectArrayElement(resultArray, 0, magniArray);
        env->SetObjectArrayElement(resultArray, 1, thetaArray);

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
        jsize length = width * height;
        jbyte *nv21Byte = env->GetByteArrayElements(nv21Array, nullptr);

        // Convert Data to Int32
        intArray nv21Int(length, 0);
        for (int i = 0; i < length; ++i) nv21Int[i] = static_cast<int>(nv21Byte[i] & 0xFF);

        // Define Gaussian Filter
        const int kernelSize = 3;
        int **gX = new int *[kernelSize];
        int **gY = new int *[kernelSize];
        const int valuesX[kernelSize][kernelSize] = {
                {-1, 0, 1},
                {-2, 0, 2},
                {-1, 0, 1}
        };
        const int valuesY[kernelSize][kernelSize] = {
                {-1, -2, -1},
                {0,  0,  0},
                {1,  2,  1}
        };
        for (int i = 0; i < kernelSize; ++i) {
            gX[i] = new int[kernelSize];
            gY[i] = new int[kernelSize];
            for (int j = 0; j < kernelSize; ++j) {
                gX[i][j] = valuesX[i][j];
                gY[i][j] = valuesY[i][j];
            }
        }

        // Perform 2D Convolution
        intArray magniInt(length, 0), thetaInt(length, 0);
        convolve2DSobel(nv21Int, magniInt, thetaInt, width, height, gX, gY, kernelSize, 4);

        // Convert Back to UInt8
        jbyteArray magniArray = env->NewByteArray(length * 3 / 2);
        jbyteArray thetaArray = env->NewByteArray(length * 3 / 2);

        jbyte *magniByte = env->GetByteArrayElements(magniArray, nullptr);
        jbyte *thetaByte = env->GetByteArrayElements(thetaArray, nullptr);
        std::memset(magniByte, 128, length * 3 / 2);
        std::memset(thetaByte, 128, length * 3 / 2);

        for (int i = 0; i < length; ++i) {
            magniByte[i] = static_cast<jbyte>(magniInt[i] & 0xFF);
            thetaByte[i] = static_cast<jbyte>(thetaInt[i] & 0xFF);
        }

        // Release Resources
        nv21Int.clear(), magniInt.clear(), thetaInt.clear();
        env->ReleaseByteArrayElements(nv21Array, nv21Byte, 0);
        env->ReleaseByteArrayElements(magniArray, magniByte, 0);
        env->ReleaseByteArrayElements(thetaArray, thetaByte, 0);


        jobjectArray resultArray = env->NewObjectArray(2,
                                   env->FindClass("[B"), nullptr);
        env->SetObjectArrayElement(resultArray, 0, magniArray);
        env->SetObjectArrayElement(resultArray, 1, thetaArray);

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
        jsize length = width * height;
        jbyte* nv21Byte = env->GetByteArrayElements(nv21Array, nullptr);

        // Convert Data to Int32
        intArray nv21Int(length, 0), tempInt(length, 0), gradInt(length, 0);
        for (int i = 0; i < length; ++i) nv21Int[i] = static_cast<int>(nv21Byte[i] & 0xFF);

        // Define Gaussian Filter
        const int kernelSize = 5;
        int* gauss_filter = new int[kernelSize];
        const int values[kernelSize] = {1, 4,  6,  4,  1};
        for (int i = 0; i < kernelSize; ++i) {
            gauss_filter[i] = values[i];
        }

        // Perform 1D Convolution
        convolve1DCol(nv21Int, tempInt, width, height, gauss_filter, kernelSize, 16);
        convolve1DRow(tempInt, gradInt, width, height, gauss_filter, kernelSize, 16);

        // Convert Back to UInt8
        jbyteArray gradArray = env->NewByteArray(length * 3/2);

        jbyte *gradByte = env->GetByteArrayElements(gradArray, nullptr);
        std::memset(gradByte, 128, length * 3/2);

        for (int i = 0; i < length; ++i) {
            gradByte[i] = static_cast<jbyte>(gradInt[i] & 0xFF);
        }

        // Release Resources
        nv21Int.clear(), tempInt.clear(), gradInt.clear();
        env->ReleaseByteArrayElements(nv21Array, nv21Byte, 0);
        env->ReleaseByteArrayElements(gradArray, gradByte, 0);

        return gradArray;

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
        jsize length = width * height;
        jbyte* nv21Byte = env->GetByteArrayElements(nv21Array, nullptr);

        // Convert Data to Int32
        intArray nv21Int(length, 0), gradInt(length, 0);
        for (int i = 0; i < length; ++i) nv21Int[i] = static_cast<int>(nv21Byte[i] & 0xFF);

        // Define Gaussian Filter
        const int kernelSize = 5;
        int** kernel = new int*[kernelSize];
        const int values[kernelSize][kernelSize] = {
                {1, 4,  6,  4,  1},
                {4, 16, 24, 16, 4},
                {6, 24, 36, 24, 6},
                {4, 16, 24, 16, 4},
                {1, 4,  6,  4,  1}
        };
        for (int i = 0; i < kernelSize; ++i) {
            kernel[i] = new int[kernelSize];
            for (int j = 0; j < kernelSize; ++j) {
                kernel[i][j] = values[i][j];
            }
        }

        // Perform 2D Convolution
        convolve2D(nv21Int, gradInt, width, height, kernel, kernelSize, 256);

        // Convert Back to UInt8
        jbyteArray gradArray = env->NewByteArray(length * 3/2);

        jbyte *gradByte = env->GetByteArrayElements(gradArray, nullptr);
        std::memset(gradByte, 128, length * 3/2);

        for (int i = 0; i < length; ++i) {
            gradByte[i] = static_cast<jbyte>(gradInt[i] & 0xFF);
        }

        // Release Resources
        nv21Int.clear(), gradInt.clear();
        env->ReleaseByteArrayElements(nv21Array, nv21Byte, 0);
        env->ReleaseByteArrayElements(gradArray, gradByte, 0);

        return gradArray;

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
        jsize length = width * height;
        auto magniArray = (jbyteArray) env->GetObjectArrayElement(sobelArray, 0);
        auto thetaArray = (jbyteArray) env->GetObjectArrayElement(sobelArray, 0);

        jbyte* magniByte = env->GetByteArrayElements(magniArray, nullptr);
        jbyte* thetaByte = env->GetByteArrayElements(thetaArray, nullptr);

        // Convert Data to Int32
        intArray magniInt(length, 0), thetaInt(length, 0), nmsInt(length, 0);
        for (int i = 0; i < length; ++i) {
            magniInt[i] = static_cast<int>(magniByte[i] & 0xFF);
            thetaInt[i] = static_cast<int>(thetaByte[i] & 0xFF);
        }

        // Perform Non-Max Suppression
        nonMaxSuppression(magniInt, thetaInt, nmsInt, width, height);

        // Convert Back to UInt8
        jbyteArray nmsArray = env->NewByteArray(length * 3/2);

        jbyte *nmsByte = env->GetByteArrayElements(nmsArray, nullptr);
        std::memset(nmsByte, 128, length * 3/2);

        for (int i = 0; i < length; ++i) {
            nmsByte[i] = static_cast<jbyte>(nmsInt[i] & 0xFF);
        }

        // Release Resources
        magniInt.clear(), thetaInt.clear(), nmsInt.clear();
        env->ReleaseByteArrayElements(nmsArray, nmsByte, 0);
        env->ReleaseByteArrayElements(magniArray, magniByte, 0);
        env->ReleaseByteArrayElements(thetaArray, thetaByte, 0);

        return nmsArray;

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
        jsize length = width * height;
        jbyte* nmsByte = env->GetByteArrayElements(nmsArray, nullptr);

        // Convert Data to Int32
        intArray nmsInt(length, 0), threshInt(length, 0);
        for (int i = 0; i < length; ++i) nmsInt[i] = static_cast<int>(nmsByte[i] & 0xFF);

        thresholding(nmsInt, threshInt, width, height, weak_pixel, strong_pixel,
                     low_thresh_ratio, high_thresh_ratio);

        // Convert Back to UInt8
        jbyteArray threshArray = env->NewByteArray(length * 3/2);

        jbyte *threshByte = env->GetByteArrayElements(threshArray, nullptr);
        std::memset(threshByte, 128, length * 3/2);

        for (int i = 0; i < length; ++i) {
            threshByte[i] = static_cast<jbyte>(threshInt[i] & 0xFF);
        }

        // Release Resources
        nmsInt.clear(), threshInt.clear();
        env->ReleaseByteArrayElements(nmsArray, nmsByte, 0);
        env->ReleaseByteArrayElements(threshArray, threshByte, 0);

        return threshArray;

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
        jsize length = width * height;
        jbyte* threshByte = env->GetByteArrayElements(threshArray, nullptr);

        // Convert Data to Int32
        intArray threshInt(length, 0), cannyInt(length, 0);
        for (int i = 0; i < length; ++i) threshInt[i] = static_cast<int>(threshByte[i] & 0xFF);

        hysteresis(threshInt, cannyInt, width, height, weak_pixel, strong_pixel);

        // Convert Back to UInt8
        jbyteArray cannyArray = env->NewByteArray(length * 3/2);

        jbyte *cannyByte = env->GetByteArrayElements(cannyArray, nullptr);
        std::memset(cannyByte, 128, length * 3/2);

        for (int i = 0; i < length; ++i) {
            cannyByte[i] = static_cast<jbyte>(cannyInt[i] & 0xFF);
        }

        // Release Resources
        threshInt.clear(), cannyInt.clear();
        env->ReleaseByteArrayElements(cannyArray, cannyByte, 0);
        env->ReleaseByteArrayElements(threshArray, threshByte, 0);

        return cannyArray;

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
        jsize length = width * height;
        jbyte* cannyByte = env->GetByteArrayElements(cannyArray, nullptr);
        jdouble *bins1 = env->GetDoubleArrayElements(bins1Array, nullptr);
        jdouble *bins2 = env->GetDoubleArrayElements(bins2Array, nullptr);

        // Convert Data to Int32
        intArray cannyInt(length, 0);
        for (int i = 0; i < length; ++i) cannyInt[i] = static_cast<int>(cannyByte[i] & 0xFF);

        floatArray houghTransformData = houghTransform(cannyInt, width, height, rho, theta, bins1, bins2, thresh);

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