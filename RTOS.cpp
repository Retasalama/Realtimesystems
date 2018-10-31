#include "mbed.h"
#include <aani.h>
#include <alipaasto_kertoimet.h>
#include <twiddle.h>

#define TS 0.000125
#define KOKO 20
#define KOKO2 16

DigitalOut myled(LED1);
AnalogIn AnalogToDigital (A0);
AnalogOut DigitalToAnalog (A2);
Serial pc(USBTX, USBRX);

typedef struct {
    int indeksi_aanitaul;
    float aanitaul[KOKO];
} message_t;

typedef struct {
    int a;
    float aaniFFTlle[KOKO2];
} mail_t;

EventFlags bitFlag32;
MemoryPool<message_t, 20> mpool;
Queue<message_t, 20> queue;
Queue<message_t, 20> queue2;
Mail<mail_t, 20> mail_box;

message_t *ptr;
message_t *writeptr;
message_t *suodatinptr;
mail_t *mailPointteri;
Ticker kello;
bool mullaOnPointteri;
bool mullaOnWritePointteri;
int indeksi;
float muisti[3][3]={{0}};
float tulos_cos[16]={0};
float tulos_sin[16]={0};
float tulos_abs[16]={0};
float input[16]={0};


void filterFunction(void);
void fftFunction(void);
float IIR_suodatin(float nayte,int suodatin);
void laskeFft(float * input, float * real, float * imag, float * abs);

void kellonpalvelija(void){
    
    if (mullaOnPointteri == false) {
        ptr = mpool.alloc();
        mullaOnPointteri = true;
        ptr->indeksi_aanitaul = 0;
    }
    
    ptr->aanitaul[ptr->indeksi_aanitaul] = taulukko[indeksi];
    ptr->indeksi_aanitaul++;
    
    if(ptr->indeksi_aanitaul == 20){
        queue.put(ptr);
        mullaOnPointteri = false;
    }
    
    indeksi = (indeksi+1)%80000;
    
    if(mullaOnWritePointteri == false){
        osEvent evt = queue2.get(0);
        if(evt.status == osEventMessage) {
            writeptr = (message_t*)evt.value.p;
            mullaOnWritePointteri = true;
            writeptr->indeksi_aanitaul = 0;
        }
    }
    if(mullaOnWritePointteri == true) {
        DigitalToAnalog = writeptr->aanitaul[writeptr->indeksi_aanitaul];
        writeptr->indeksi_aanitaul++;
        
        if(writeptr->indeksi_aanitaul == 20) {
            mullaOnWritePointteri = false;
            mpool.free(writeptr);
        }
    }
    
}

int main(){
    indeksi = 0;
    mullaOnPointteri = false;
    mullaOnWritePointteri = false;
    kello.attach(&kellonpalvelija, TS);
    Thread filterThread(osPriorityNormal);
    Thread fftThread(osPriorityNormal);
    
    filterThread.start(filterFunction);
    fftThread.start(fftFunction);
    while(true){
        //myled = 1;
        //wait(0.5);
        //myled = 0;
        //pc.printf("jep");
    }
}

void filterFunction(void) {
    while(true){
        osEvent evt = queue.get();
        if(evt.status == osEventMessage) {
            suodatinptr = (message_t*)evt.value.p;
            uint32_t lipun_arvo = bitFlag32.get();
            if(lipun_arvo&0x1) {
                mail_t *mail = mail_box.alloc();
                for(int j = 0; j<KOKO2; j++) {
                    mail->aaniFFTlle[j] = suodatinptr->aanitaul[j];
                }
                mail_box.put(mail);
                bitFlag32.clear();
            }
            
            for(int i = 0; i>KOKO; i++) {
                float nayte = suodatinptr->aanitaul[i];
                nayte = IIR_suodatin(nayte,0);
                nayte = IIR_suodatin(nayte,1);
                nayte = IIR_suodatin(nayte,2);
                suodatinptr->aanitaul[i] = nayte;
            }
            queue2.put(suodatinptr);
        }
    }
}

void fftFunction(void) {
    while(true){
        bitFlag32.set(0x1);
        osEvent evt = mail_box.get();
        if(evt.status == osEventMail) {
            mailPointteri = (mail_t*)evt.value.p;
            for(int i = 0; i<KOKO2; i++) {
                //suodattamattomat näytteet input-taulukkoon
                input[i] = mailPointteri->aaniFFTlle[i];
            }
            //fft:n laskeminen
            laskeFft(input, tulos_cos, tulos_sin, tulos_abs);
            //real, imag ja abs-tulostus
            for(int j = 0 ; j<KOKO2; j++) {
                printf("Real: %f\n\r", tulos_cos[j]);
                printf("Imag: %f\n\r", tulos_sin[j]);
                printf("Abs: %f\n\n\r", tulos_abs[j]);
            }
            mail_box.free(mailPointteri);
        }
    }
}


float IIR_suodatin(float nayte,int suodatin)
{
    float tulos = 0;
    
    muisti[suodatin][2]=muisti[suodatin][1];
    muisti[suodatin][1]=muisti[suodatin][0];
    muisti[suodatin][0]= nayte - SOS[suodatin][4]*muisti[suodatin][1] - SOS[suodatin][5]*muisti[suodatin][2];
    
    tulos = SOS[suodatin][0]*muisti[suodatin][0] + SOS[suodatin][1]*muisti[suodatin][1] + SOS[suodatin][2]*muisti[suodatin][2];
    
    return tulos;
}

void laskeFft(float * input, float * real, float * imag, float * abs)
{
    int N = 16;
    for(int I = 0; I < N ; I++)
    {
        real[I]= 0;
        imag[I]= 0;
    }
    
    //for(int m = 1 ; m<6 ; m++) // lasketaan vain 1 = 500, 2=1000, 3=1500, 4=2000, 5 = 2500 Hz
    for(int m = 0 ; m<16 ; m++) // lasketaan kaikki taajuuskomponentit
    {
        for(int n = 0 ; n < N ; n++)
        {
            real[m] = real[m] + input[n]*kosini[m][n];
            imag[m] = imag[m] + input[n]*sini[m][n];
            
        }
        abs[m] = sqrt(real[m]*real[m] + imag[m]*imag[m]);
    }
}



