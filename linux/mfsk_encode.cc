#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mfsk.h"

static char *InputFileName=0;        // the name of the input file (specified on the command line)
static char *AudioFileName;          // file name to save audio
static FILE *InputFile=0;
static FILE *AudioFile=0;

static MFSK_Parameters<float>  Parameters;
static MFSK_Transmitter<float> Transmitter;

int main(int argc, char *argv[])
{ int Error;

  // Parameters.Default();
  // Parameters.Preset();

  // read the command line options and the input file name (if specified)
  int arg;
  int Help=0;
  for(arg=1; arg<argc; arg++)
  { if(argv[arg][0]=='-') // if '-' then this is an option, otherwise the name of a file
    { int Error=Parameters.ReadOption(argv[arg]);
      if(Error<0)
      { printf("Invalid parameter(s) in %s\n",argv[arg]); }
      else if(Error==0)
      { Help=1; }
    }
    else // if does not start with a dash, then it must be a file name
    {      if(InputFileName==0) InputFileName=argv[arg]; // assign an text input
      else if(AudioFileName==0) AudioFileName=argv[arg]; // assign as audio output
      else Help=1;
    }
  }

  if(Help)
  { printf("\nmfsk_encode [options] <text file> <audio file>\nPptions:\n");
    printf("%s\n",Parameters.OptionHelp());
    return -1; }

  if(InputFileName)
  { InputFile=fopen(InputFileName,"rt");
    if(InputFile==0)
    { printf("Can not open %s for read !\n", InputFileName);
      return -1; }
  }

  if(AudioFileName)
  { AudioFile=fopen(AudioFileName,"wb");
    if(AudioFile==0)
    { printf("Can not open %s for write !\n", AudioFileName);
      return -1; }
  }

  if(InputFileName==0 || AudioFileName==0)
  { printf("Both input and audio file names need to be given\n"); return -1; }

  Error=Parameters.Preset();
  if(Error<0)
  { printf("Parameters.Preset() => %d\n",Error); return -1; }

  // preset the transmitter's internal arrays (negative return means fatal error)
  Error=Transmitter.Preset(&Parameters);
  if(Error<0)
  { printf("Transmitter.Preset() => %d\n",Error); return -1; }

  int16_t AudioBuffer[Transmitter.MaxOutputLen];

  printf("MFSK encoder by Pawel Jalocha, October 2023\n");
  Parameters.Print(1);

  printf("MFSK encodng %s into audio in %s\n", InputFileName, AudioFileName);

  int Head = 8000;        // [samples]
  int Tail = 8000;        // [samples]
  Transmitter.Start();
  // int Stop=0;
  for( ; ; )
  { if(!Transmitter.Running() && Head<=0 && Tail<=0) break;

    if(Head<=0)
    { int Char=fgetc(InputFile);
      if(Char!=EOF)
      { // printf("Tx: %02X\n", Char);
        Transmitter.PutChar((uint8_t)Char); }
      else Transmitter.Stop();
    }

    uint8_t TxChar;
    while(Transmitter.GetChar(TxChar)>0)          // monitor the characters being transmitted
    { printf("%c", TxChar); }                     // and write them to the screen
    // fflush(stdout);

    int Len=Transmitter.Output(AudioBuffer);      // read audio from the transmitter
    if(Len>0)
    { if(Head>0) Head-=Len;
      else if(!Transmitter.Running() && Tail>0) Tail-=Len;
      // printf("Tx: %d samples ( Head=%d, Tail=%d)\n", Len, Head, Tail);
      if(fwrite(AudioBuffer, 2, Len, AudioFile)!=Len) Transmitter.Stop(); // write the audio into the file
    }
  }

  printf("\n");

  if(AudioFile) fclose(AudioFile);
  if(InputFile) fclose(InputFile);

  return 0; }

