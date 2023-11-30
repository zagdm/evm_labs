#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

using namespace cv;

int main() {
	VideoCapture capture(0);
	if (!capture.isOpened()) {
		std::cout << "Ошибка при открытии камеры" << std::endl;
		return 1;
	}
	Mat frame;
	double fps = 0;
	double t_fps = 0;
	double t_input = 0;
	double t_proc = 0;
	double t_out = 0;
	double average = 0;
	double min = 100000;
	while (true) {
		namedWindow("Camera", WINDOW_NORMAL);
		resizeWindow("Camera", 1280, 720);
		t_fps = (double)getTickCount() / getTickFrequency();
		capture.read(frame);
		t_input = (double)getTickCount() / getTickFrequency();
		if (frame.empty()) {
			std::cout << "Ошибка чтения кадра" << std::endl;
			break;
		}
		circle(frame, Point(frame.cols / 2, frame.rows / 2), 50, Scalar(0, 255, 255), FILLED);
		rectangle(frame, Point(frame.cols / 4, frame.rows / 4), Point(3 * frame.cols / 4, 3 * frame.rows / 4), Scalar(255, 255, 0), 10);
		t_proc = (double)getTickCount() / getTickFrequency();
		putText(frame, "FPS: " + std::to_string(fps), Point(10, 20), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0, 0, 255), 2);
		putText(frame, "Input time: " + std::to_string(t_input - t_fps) + "s", Point(10, 40), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0, 0, 255), 2);
		putText(frame, "Proccesing time: " + std::to_string(t_proc - t_input) + "s", Point(250, 20), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0, 0, 255), 2);
		putText(frame, "Output time: " + std::to_string(t_out) + "s", Point(250, 40), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0, 0, 255), 2);
		t_out = (double)getTickCount() / getTickFrequency();
		imshow("Camera", frame);
		t_out = (double)getTickCount() / getTickFrequency() - t_out;
		t_fps = (double)getTickCount() / getTickFrequency() - t_fps;
		if (average) {
			average = (average + t_fps) / 2;
		} else {
			average = t_fps;
		}
		if (t_fps < min) {
			min = t_fps;
		}
		fps = 1 / t_fps;
		if (waitKey(30) == 27) {
			break;
		}
	}
	destroyWindow("Camera");
	capture.release();
	std::cout << "Минимальное время обработки: " << min << "с" << std::endl;
	std::cout << "Среднее время обработки: " << average << "с" << std::endl;
	return 0;
}
