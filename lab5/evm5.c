#include <opencv2/highgui/highgui.hpp>

int main(int argc,char *argv[]) {
	CvCapture *capture = cvCreateCameraCapture(0);
	if (!capture) return 0;
	int fps = 0; 
		double t = 0;
	while(1) {
		IplImage *frame = cvQueryFrame(capture);
		if(!frame) break;
	GaussianBlur(frame, frame, Size(5, 5), 1.5);
	t = (double)getTickCount();
		cvShowImage("test", frame);
		char c = cvWaitKey(33);
		if(c == 27) break;
		t = ((double)getTickCount() - t) / getTickFrequency();
	fps = 1 / t;
	putText(frame, "FPS: " + to_string(fps), Point(10, 30), 						FONT_HERSHEY_SIMPLEX, 1, Scalar(0, 0, 255), 2);
	}
	cvReleaseCapture(&capture);
	cvDestroyWindow("test");
}