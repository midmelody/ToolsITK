#include <stdio.h>

void main() {
	unsigned short i, pi;
	i = 0; pi = 0;
	while(true) {
		i= i + 10;
		if(i < pi) {
			printf("max value %u\n", pi);
			break;
		}
		pi = i;
	}
}
