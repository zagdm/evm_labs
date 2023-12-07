#include <libusb-1.0/libusb.h>
#include <iostream>

void print_title() {
	std::cout << "===============================================" << std::endl;
	std::cout << "| * ����� ����������" << std::endl;
	std::cout << "|   | * ����� ������ ����������" << std::endl;
	std::cout << "|   |    | * ������������� �������������" << std::endl;
	std::cout << "|   |    |     | * ������������� ����������" << std::endl;
	std::cout << "|   |    |     |    |  * �������� �����" << std::endl;
	std::cout << "+---+----+-----+----+--------------------------" << std::endl;
}

void print_devices(size_t dev_numb, libusb_device* dev) {
	printf(" %3.ld:", dev_numb);
	libusb_device_handle* handle = nullptr; 
	libusb_device_descriptor desc{};
	int x = libusb_get_device_descriptor(dev, &desc);
	if (x < 0) {
		printf("���������� ���������� �� �������. ���: %d.\n", x);
		return;
	}
	printf("%.4x %.5x %.4x ", desc.bDeviceClass, desc.idVendor, desc.idProduct);
	libusb_open(dev, &handle);
	if (handle && desc.iSerialNumber) {
		unsigned char serial_numb[256];
		x = libusb_get_string_descriptor_ascii(handle, desc.iSerialNumber, serial_numb, sizeof(serial_numb));
		std::cout << serial_numb << std::endl;
	}  else {
		std::cout << "empty" << std::endl;
	}
	libusb_close(handle);
}

int main() {
	libusb_context* context = nullptr;
	libusb_device** devices;
	int x = libusb_init(&context);
	if (x < 0) {
		fprintf(stderr, "������ �������������. ���: %d.\n", x);
		return 1;
	}
	size_t number = libusb_get_device_list(context, &devices);
	if (number < 0) {
		fprintf(stderr, "C����� ��������� �� �������. ���: %d\n", x);
		libusb_exit(context);
		return 1;
	}
	std::cout << "������� " << number << "���������" << std::endl;
	print_title();
	for (size_t i = 0; i < number; i++) print_devices(i + 1, devices[i]);
	std::cout << "===============================================" << std::endl;
	libusb_free_device_list(devices, 1);
	libusb_exit(context);
	return 0;
}