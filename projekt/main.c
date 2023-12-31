#include <stdint.h>   /* Declarations of uint_32 and the like */
#include <pic32mx.h>  /* Declarations of system-specific addresses etc */
#include "header.h"  /* Declatations for these labs */

int main(void) {

	SYSKEY = 0xAA996655;  /* Unlock OSCCON, step 1 */
	SYSKEY = 0x556699AA;  /* Unlock OSCCON, step 2 */
	while(OSCCON & (1 << 21)); /* Wait until PBDIV ready */
	OSCCONCLR = 0x180000; /* clear PBDIV bit <0,1> */
	while(OSCCON & (1 << 21));  /* Wait until PBDIV ready */
	SYSKEY = 0x0;  /* Lock OSCCON */
	
	/* Set up output pins */
	AD1PCFG = 0xFFFF;
	ODCE = 0x0;
	TRISECLR = 0xFF;
	PORTE = 0x0;
	
	/* Output pins for display signals */
	PORTF = 0xFFFF;
	PORTG = (1 << 9);
	ODCF = 0x0;
	ODCG = 0x0;
	TRISFCLR = 0x70;
	TRISGCLR = 0x200;
	
	/* Set up input pins */
	TRISDSET = (1 << 8);
	TRISFSET = (1 << 1);
	
	/* Set up SPI as master */
	SPI2CON = 0;
	SPI2BRG = 4;
	/* SPI2STAT bit SPIROV = 0; */
	SPI2STATCLR = 0x40;
	/* SPI2CON bit CKP = 1; */
        SPI2CONSET = 0x40;
	/* SPI2CON bit MSTEN = 1; */
	SPI2CONSET = 0x20;
	/* SPI2CON bit ON = 1; */
	SPI2CONSET = 0x8000;
	
	display_init();

	display_update();

    float time=0.0f;
    inputInit();
		int choice=0;
		int choiceMade= 1;
    while(choiceMade){
		volatile unsigned int pointerToPORTE;
		pointerToPORTE = (volatile unsigned int)0xbf886110;
		unsigned int currBtns = getbtns();

		display_update();
        display_string( 0, "Btn1 for cube" );
        display_string( 1, "Btn2 for pyramid" );
        display_string( 2, "Btn3 for block" );
		display_string( 3, "Rst by btn Reset" );
        if (currBtns) {
            if (currBtns & 0x1)  // If bit number one is set, decrement fSigma
            {
                choice=0;
				choiceMade=0;
            }
            if (currBtns & 0x2)  // If bit number two is set, increment fTheta
            {
                choice=1;
				choiceMade=0;
            }
            if (currBtns & 0x4)  // If bit number four is set, increment fSigma
            {
                choice=2;
				choiceMade=0;
            }
			if (currBtns & 0x8)  // If bit number four is set, increment fSigma
            {
                choice=3;
				choiceMade=0;
            }
    }
	}
	display_string( 0, "" );
	display_string( 1, "Loading..." );
	display_string( 2, "" );
	if (choice==3){
		display_string( 2, "secret 4th" );
	}
	
	display_string( 3, "" );
	display_update();
    labinit(choice);
	display_string( 2, "" );
	display_string( 1, "" );
    while(1)
    {
          newlabwork(); /* Do lab-specific things again and again */
    }
    return 0;

}

