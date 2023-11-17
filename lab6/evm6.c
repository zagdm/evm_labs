#include <libusb-1.0/libusb.h>
#include <iostream>

void print_title() {
    printf("===================================\n");
    printf("|* класс устройства\n");
    printf("|  | * идентификатор производителя\n");
    printf("|  |    | * идентификатор устройства\n");
    printf("|  |    |    | * серийный номер\n");
    printf("+--+----+----+--------------------+\n");
}

void print_devices(libusb_device* dev) {
    libusb_device_handle* handle = NULL; 
    libusb_device_descriptor desc{};
    int x = libusb_get_device_descriptor(dev, &desc);
    if (x < 0) {
        fprintf(stderr, "Дескриптор устройства не получен. Код: %d.\n", x);
        return;
    }
    printf("%.2x %.4x %.4x ", desc.bDeviceClass, desc.idVendor, desc.idProduct);
    libusb_open(dev, &handle);
    if (handle && desc.iSerialNumber) {
        char serial_numb[256];
        x = libusb_get_string_descriptor_ascii(handle, desc.iSerialNumber, serial_numb, sizeof(serial_numb));
        printf("%s\n", serial_numb);
    }  else {
        printf("empty\n");
    }
}

int main() {
    libusb_device** devices;
    libusb_context* context = NULL;
    size_t number;
    int x = libusb_init(&context);
    if (x < 0) {
        fprintf(stderr, "Ошибка инициализации. Код: %d.\n", x);
        return 1;
    }
    size_t number = libusb_get_device_list(context, &devices);
    if (number < 0) {
        fprintf(stderr, "Cписок устройств не получен. Код: %d\n", x);
        return 1;
    }
    printf("Найдено %ld устройств\n", number);
    print_title();
    for (size_t i = 0; i < number; i++) {
        print_devices(devices[i]);
    }
    printf("===================================\n");
    libusb_free_device_list(devices, 1);
    libusb_exit(context);
    return 0;
}