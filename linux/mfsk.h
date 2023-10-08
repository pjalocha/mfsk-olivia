// MFSK trasnmitter and receiver code,
// (c) Pawel Jalocha, April 2006

#ifndef __MFSK_H__
#define __MFSK_H__

// =====================================================================

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "struc.h"
#include "fht.h"
#include "cmpx.h"
#include "fft.h"
#include "gray.h"
#include "lowpass3.h"
#include "buffer.h"
#include "rateconv.h"

#include "noise.h"

// =====================================================================
/*

The convention to pass parameters to the objects:

User-setable parameters are listed on top of the class.
When an object is being created it is given certain default parameters
by the Default() call. If the user wishes to modify some (or all) of them
he should refer to them directly. Then, the user must call Preset()
which will preset the internal object structures including dependend parameters
and dynamic arrays. Only then calls line Input() and Process() can be executed.

Preset() may return a negative number in case some parameter are not valid.
It can as well adjust some parameters to the closest valid values.

If the user wishes to change the parameters he should do so and then call Preset() again,
however, the data accumulated in the internal structures are lost.

When the user wants to save memory and free the internal storage allocated
by an object he can call Free() and then later he can call Preset() again to re-use
this object.

*/
// =====================================================================

// fast (integer) power of two
static inline int Exp2(uint32_t X)
{ return (uint32_t)1<<X; }

// fast (integer) base-2 logarythm
static inline int Log2(uint32_t X)
{ uint32_t Y;
  for( Y=0; X>1; X>>=1)
    Y+=1;
  return Y; }

// =====================================================================

template <class Type>
 int FitPeak(Type &PeakPos, Type &Peak, Type Left, Type Center, Type Right)
{ Type A=(Right+Left)/2-Center;
  if(A>=0) return -1;
  Type B=(Right-Left)/2;
  PeakPos=(-B/(2*A));
  Peak = A*PeakPos*PeakPos + B*PeakPos + Center;
  return 1; }

// =====================================================================

template <class Type>
 Type Limit(Type X, Type Lower, Type Upper)
{ if(X>Upper) return Upper;
  if(X<Lower) return Lower;
  return X; }

// =====================================================================

// convert audio from floating point to 16-bit signed
template <class Type>
 void ConvertToS16(Type *Input, int16_t *Output, int Len, Type Scale=32768.0)
{ int Idx;
  const int Limit=32767;
  for(Idx=0; Idx<Len; Idx++)
  { int Out=(int)floor(Scale*Input[Idx]+0.5);
    if(Out>Limit) Out=Limit;
	else if(Out<(-Limit)) Out=(-Limit);
    Output[Idx]=Out; }
}

template <class Type>
 int ConvertToS16(Seq<Type> &Input, Seq<int16_t> &Output, Type Scale=32768.0)
{ if(Output.EnsureSpace(Input.Len)<0) return -1;
  Output.Len=Input.Len;
  ConvertToS16(Input.Elem, Output.Elem, Input.Len, Scale);
  return 0; }

// =====================================================================

// =====================================================================

// the symbol shape described in frequency domain
static const double MFSK_SymbolFreqShape[] =
  { +1.0000000000, +2.1373197349, +1.1207588117, -0.0165609232 } ;
static const int MFSK_SymbolFreqShapeLen=sizeof(MFSK_SymbolFreqShape)/sizeof(double);

// the basic parameters
template <class FloatType=float>
 class MFSK_Parameters
{ public:
                                             // primary parameters
        int BitsPerSymbol;                   // [Bits] number of bits per MFSK symbol/tone
        int Bandwidth;                       // [Hz]   total bandwith like 500, 1000 or 2000Hz
        int SampleRate;                      // [Hz]   audio sample rate: 8000Hz by defaults
  FloatType LowerBandEdge;                   // [Hz]   
  FloatType InputSampleRate;                 // [Hz]   true audio card input sample rate, can be different from the ideal 8kHz
  FloatType OutputSampleRate;                // [Hz]   true audio card output sample rate
                                             // receiver parameters
        int RxSyncMargin;                    // [MFSK carriers] lock search frequency margin
        int RxSyncIntegLen;                  // [FEC Blocks]    averaging size for the lock search
  FloatType RxSyncThreshold;                 // [S/N]           minimal SNR for a lock and decoding

                                             // fixed parameters
  static const int BitsPerCharacter   = 7;   // [Bits]
  static const int SymbolsPerBlock    = 1<<(BitsPerCharacter-1);
  static const int CarrierSepar       = 4;   // [FFT bins]
  static const int SpectraPerSymbol   = 4;   // [Spectral (FFT) slices]
  static const int SpectraPerBlock    = SpectraPerSymbol*SymbolsPerBlock;
  static const int UseGrayCode        = 1;
  static const int PhaseDiffer        = 1;
  static const int RxSyncSquareEnergy = 1;
  static const int DecodeSquareEnergy = 1;
  static const uint64_t ScramblingCode   = 0xE257E6D0291574ECLL;

                                          // secondary parameters
  int Carriers;
  int SymbolSepar;                        // [Samples]   distance between succesive symbols
  int SymbolLen;                          // [Samples]   size of the FFT slice produce the symbol shape
  int FirstCarrier;                       // [FFT bins]

  MFSK_Parameters()
    { Default(); }

  void Default(void)
    { BitsPerSymbol       = 5;
      SampleRate          = 8000;
      Bandwidth           = 1000;
      LowerBandEdge       = SampleRate/16;
      InputSampleRate     = SampleRate;
      OutputSampleRate    = SampleRate;
      RxSyncIntegLen      = 8;
      RxSyncMargin        = 4;
      RxSyncThreshold     = 3.0; }

  int Preset(void)
    {
      if(BitsPerSymbol>8) BitsPerSymbol=8;
      else if(BitsPerSymbol<1) BitsPerSymbol=1;
      Carriers=Exp2(BitsPerSymbol);

      int MinBandwidth=SampleRate/64;
      int MaxBandwidth=SampleRate/4;
      if(Bandwidth<MinBandwidth) Bandwidth=MinBandwidth;
      else if(Bandwidth>MaxBandwidth) Bandwidth=MaxBandwidth;
      Bandwidth=MinBandwidth*Exp2(Log2(Bandwidth/MinBandwidth));

      SymbolSepar=(SampleRate/Bandwidth)*Carriers;
      SymbolLen=SymbolSepar*CarrierSepar;

      // printf("Preset() SampleRate:%d, CarrierSepar:%d, SymbolLen:%d, FirstCarrier:%d, LowerBandEdge:%3.1f Hz\n",
      //                   SampleRate, CarrierSepar, SymbolLen, FirstCarrier, LowerBandEdge);

      FirstCarrier=(int)floor((LowerBandEdge/SampleRate)*SymbolLen+0.5)+(CarrierSepar/2);
      if((FirstCarrier+Carriers*CarrierSepar)>=(SymbolLen/2))
        FirstCarrier=(SymbolLen/2)-Carriers*CarrierSepar;
      LowerBandEdge=(FloatType)(FirstCarrier-CarrierSepar/2)*SampleRate/SymbolLen;

      if(RxSyncMargin>(FirstCarrier/CarrierSepar)) RxSyncMargin=(FirstCarrier/CarrierSepar);

      // printf("Preset() SampleRate:%d, CarrierSepar:%d, SymbolLen:%d, FirstCarrier:%d, LowerBandEdge:%3.1f Hz\n",
      //                   SampleRate, CarrierSepar, SymbolLen, FirstCarrier, LowerBandEdge);

      return 0; }

   const char *OptionHelp(void)
     { return "\
  -T<tones>             number of tones: 4, 8, 16, [32], 64, 128, 256\n\
  -B<bandwidth>/<edge>  bandwidth: 125, 250, 500, [1000], 2000\n\
                        and lower audio band edge [500] [Hz]\n\
  -S<threshold>         S/N threshold [3.0]\n\
  -M<margin>            frequency search margin [4]\n\
  -I<period>            synchr. integration period [8]\n\
  -R<Tx>/<Rx>           the true sample rate for Tx and Rx [8000.0/8000.0]\n\
";   }

   int ReadOption(const char *Option)
     { if(Option[0]!='-') return 0;
       switch(Option[1])
       { case 'T':
          int Tones;
          if(sscanf(Option+2,"%d",&Tones)==1)
          { BitsPerSymbol=Log2(Tones); }
          else return -1;
          break;
         case 'B':
          int Band; float Edge;
          if(sscanf(Option+2,"%d/%f",&Band,&Edge)==2)
          { Bandwidth=Band; LowerBandEdge=Edge; }
          else if(sscanf(Option+2,"%d",&Band)==1)
          { Bandwidth=Band; }
          else return -1;
          break;
         case 'M':
          int Margin;
          if(sscanf(Option+2,"%d",&Margin)==1)
          { RxSyncMargin=Margin; }
          else return -1;
          break;
         case 'I':
          int IntegLen;
          if(sscanf(Option+2,"%d",&IntegLen)==1)
          { RxSyncIntegLen=IntegLen; }
          else return -1;
          break;
         case 'S':
          float Threshold;
          if(sscanf(Option+2,"%f",&Threshold)==1)
          { RxSyncThreshold=Threshold; }
          else return -1;
          break;
         case 'R':
          float SampleRate_Out,SampleRate_Inp;
          if(sscanf(Option+2,"%f/%f",&SampleRate_Out,&SampleRate_Inp)==2)
          { InputSampleRate=SampleRate_Inp; OutputSampleRate=SampleRate_Out; }
          else if(sscanf(Option+2,"%f",&SampleRate_Out)==1)
          { InputSampleRate=SampleRate_Out; InputSampleRate=SampleRate_Out; }
		  else return -1;
          break;
	     default: return 0; }
       return 1; }

   void Print(bool TxOnly=0) const
   { printf("MFSK_Parameters:\n");
     printf("%d (%4.1f-%4.1f) Hz, %d tones\n",
	           Bandwidth, LowerBandEdge, LowerBandEdge+Bandwidth, Carriers);
     printf("Sample rate: %d(int.) %6.1f(input) %6.1f(output) [Hz]\n",
	           SampleRate, InputSampleRate, OutputSampleRate);
     printf("Symbol/FFT: %d/%d, FirstCarrier=%d, FFT sampling [TxF]: %dx%d\n",
	           SymbolSepar, SymbolLen, FirstCarrier, SpectraPerSymbol, CarrierSepar);
     printf("%d bits/symbol, %5.3f baud, %d symbols/block, %5.3f sec/block\n",
	           BitsPerSymbol, BaudRate(), SymbolsPerBlock, BlockPeriod() );
     if(!TxOnly)
       printf("Synchronizer: +/-%d carriers = +/-%4.1f Hz,  %d blocks = %3.1f sec\n",
	           RxSyncMargin, RxSyncMargin*CarrierBandwidth(), RxSyncIntegLen, RxSyncIntegLen*BlockPeriod() );
   }

   FloatType BaudRate(void) const
   { return (FloatType)SampleRate/SymbolSepar; }

   FloatType FFTbinBandwidth(void) const
   { return (FloatType)SampleRate/SymbolLen; }

   FloatType CarrierBandwidth(void) const
   { return (FloatType)SampleRate/SymbolLen*CarrierSepar; }

   FloatType TuneMargin(void) const
   { return CarrierBandwidth()*RxSyncMargin; }

   FloatType BlockPeriod(void) const
   { return (SymbolsPerBlock*SymbolSepar)/(FloatType)SampleRate; }

   FloatType CharactersPerSecond(void) const
   { return BitsPerSymbol*(FloatType)SampleRate/(SymbolsPerBlock*SymbolSepar); }

} ;

// =====================================================================
// Soft-demodulate an MFSK symbol

template <class EnergyType, class SymbolType>
 void MFSK_SoftDemodulate(SymbolType *Symbol, EnergyType *SpectraEnergy,
                          int BitsPerSymbol, int CarrierSepar=1, int UseGrayCode=1, int SquareEnergy=1)
{ int Bit,Idx;

  for(Bit=0; Bit<BitsPerSymbol; Bit++)
    Symbol[Bit]=0;
  int Carriers=Exp2(BitsPerSymbol);

  EnergyType TotalEnergy=0;
  int Freq=0;
  for(Idx=0; Idx<Carriers; Idx++)                  // loop over carriers
  { uint8_t SymbIdx=Idx;
    if(UseGrayCode) SymbIdx=BinaryCode(SymbIdx);
    EnergyType Energy=SpectraEnergy[Freq];          // energy for given carrier
    if(SquareEnergy) Energy*=Energy;               // square the energy (works better, but why ?)
    TotalEnergy+=Energy;
    uint8_t Mask=1;
    for(Bit=0; Bit<BitsPerSymbol; Bit++)           // soft decision (contribution) for every bit
    { if(SymbIdx&Mask) Symbol[Bit]-=Energy;        // add or subtract the contribution
                  else Symbol[Bit]+=Energy;        // depending on bit value
      Mask<<=1; }
    Freq+=CarrierSepar; }

  if(TotalEnergy>0)                                 // normalize the soft decisions
  { for(Bit=0; Bit<BitsPerSymbol; Bit++)
      Symbol[Bit]/=TotalEnergy;
  }
}


template <class EnergyType, class SymbolType>
 void MFSK_SoftModulate(EnergyType *CarrierProb, SymbolType *Symbol,
                         int BitsPerSymbol, int UseGrayCode=1)
{
  int Carriers=Exp2(BitsPerSymbol);

  int Idx;
  for(Idx=0; Idx<Carriers; Idx++)                  // loop over carriers
  { uint8_t SymbIdx=Idx;
    if(UseGrayCode) SymbIdx=BinaryCode(SymbIdx);
    EnergyType Prob=1;
    int Bit;
    uint8_t Mask=1;
    for(Bit=0; Bit<BitsPerSymbol; Bit++)
	{ EnergyType BitProb=1.0;
	  if(SymbIdx&Mask) BitProb-=Symbol[Bit];
		          else BitProb+=Symbol[Bit];
      Prob*=(BitProb/2);
	  Mask<<=1; }
	CarrierProb[Idx]=Prob;
  }
}

// =====================================================================
// MFSK modulator, synthesis of the MFSK signal

template <class Type=float>
 class MFSK_Modulator
{ public:

   MFSK_Parameters<Type> *Parameters;

  public:

   int OutputLen;      // output length per transmitted symbol [samples]

  private:

   int SymbolLen;
   int SymbolSepar;

   Type *CosineTable;     // Cosine table for fast cos/sin calculation
   Type *SymbolShape;     // the shape of the symbol
   int   SymbolPhase;     // the phase of the tone being transmitted
   Type  *OutTap;         // output tap (buffer)
   int TapPtr;
   int WrapMask;

  public:

   MFSK_Modulator()
     { Init(); }

   ~MFSK_Modulator()
     { Free(); }

   void Init(void)
     { CosineTable=0;
	   SymbolShape=0;
       OutTap=0; }

   void Free(void)
     { free(CosineTable); CosineTable=0;
       free(SymbolShape); SymbolShape=0;
       free(OutTap); OutTap=0; }

   int Preset(MFSK_Parameters<Type> *NewParameters)
     { Parameters=NewParameters;

       SymbolLen=Parameters->SymbolLen;
       SymbolSepar=Parameters->SymbolSepar;

	   int Idx;

       if(ReallocArray(&CosineTable,SymbolLen)<0) goto Error;
       for(Idx=0; Idx<SymbolLen; Idx++)
         CosineTable[Idx]=cos((2*M_PI*Idx)/SymbolLen);

       if(ReallocArray(&SymbolShape,SymbolLen)<0) goto Error;

       { int Time;
         double Ampl=MFSK_SymbolFreqShape[0];
         for(Time=0; Time<SymbolLen; Time++)
	       SymbolShape[Time]=Ampl;
	   }
       int Freq;
	   for(Freq=1; Freq<MFSK_SymbolFreqShapeLen; Freq++)
       { int Time;
         double Ampl=MFSK_SymbolFreqShape[Freq];
         if(Freq&1) Ampl=(-Ampl);
		 int Phase=0;
         for(Time=0; Time<SymbolLen; Time++)
	     { SymbolShape[Time]+=Ampl*CosineTable[Phase];
           Phase+=Freq; if(Phase>=SymbolLen) Phase-=SymbolLen; }
       }
       { int Time;
         double Scale=1.0/(2*Parameters->CarrierSepar);
         for(Time=0; Time<SymbolLen; Time++)
	       SymbolShape[Time]*=Scale;
	   }

       if(ReallocArray(&OutTap,SymbolLen)<0) goto Error;
       for(Idx=0; Idx<SymbolLen; Idx++)
         OutTap[Idx]=0;
       TapPtr=0;

       WrapMask=SymbolLen-1;
       SymbolPhase=0;

       OutputLen=SymbolSepar;

       return 0;

       Error: Free(); return -1; }

   void Send(uint8_t Symbol)
     {
       if(Parameters->UseGrayCode) Symbol=GrayCode(Symbol);

       int SymbolFreq=Parameters->FirstCarrier+Parameters->CarrierSepar*Symbol;

       int TimeShift=SymbolSepar/2-SymbolLen/2;
       SymbolPhase+=SymbolFreq*TimeShift;
       SymbolPhase&=WrapMask;

       AddSymbol(SymbolFreq,SymbolPhase);

       TimeShift=SymbolSepar/2+SymbolLen/2;
       SymbolPhase+=SymbolFreq*TimeShift;
       SymbolPhase&=WrapMask;

       if(Parameters->PhaseDiffer)
       { int PhaseShift=SymbolLen/4;
         if(rand()&1) PhaseShift=(-PhaseShift);
         SymbolPhase+=PhaseShift; }

       SymbolPhase&=WrapMask;

     }

   // get output as 16-bit signed data
   int Output(int16_t *Buffer)
     { const Type Scale=32768.0;
       const int32_t Limit=0x7FFF;
       int Idx;

       for(Idx=0; Idx<SymbolSepar; Idx++)
       { Type Ampl=OutTap[TapPtr];
         Ampl*=Scale;
         int32_t Out=(int32_t)floor(Ampl+0.5);
         if(Out>Limit) Out=Limit;
         else if(Out<(-Limit)) Out=(-Limit);
         Buffer[Idx]=(int16_t)Out;
         OutTap[TapPtr]=0;
         TapPtr+=1; TapPtr&=WrapMask; }

       return SymbolSepar; }

   template <class OutType>
    int Output(OutType *Buffer)
     {  int Idx;
       for(Idx=0; Idx<SymbolSepar; Idx++)
       { Buffer[Idx]=OutTap[TapPtr];
         OutTap[TapPtr]=0;
         TapPtr+=1; TapPtr&=WrapMask; }
       return SymbolSepar; }

  private:

   void AddSymbol(int Freq, int Phase)
     { int Time;
       for(Time=0; Time<SymbolLen; Time++)
       { OutTap[TapPtr]+=CosineTable[Phase]*SymbolShape[Time];
         Phase+=Freq; Phase&=WrapMask;
         TapPtr+=1; TapPtr&=WrapMask; }
     }

} ;

// =====================================================================

// A running-box, low pass filter
template <class TapType=float, class OutType=double>
 class BoxFilter
{ public:

   int Len;

   TapType *Tap;
   int Ptr;
   OutType Output;

  BoxFilter()
    { Tap=0; }

  ~BoxFilter()
    { free(Tap); }

  void Free(void)
    { free(Tap); Tap=0; }

  int Preset(void)
    { if(ReallocArray(&Tap,Len)<0) return -1;
      Clear();
      return 0; }

  void Clear(void)
    { int Idx;
      for(Idx=0; Idx<Len; Idx++)
        Tap[Idx]=0;
      Ptr=0;
      Output=0; }

  template <class InpType>
   void Process(InpType Input)
    { Output-=Tap[Ptr];
      Output+=Input;
      Tap[Ptr]=Input;
      Ptr+=1; if(Ptr>=Len) Ptr-=Len; }

} ;

// =====================================================================

// Input processor, removes coherent interference and pulse noise
template <class Type>
 class MFSK_InputProcessor
{ public:

   // the user-settable parameters:
   int WindowLen;  // spectral analysis (FFT) window length
   Type LimiterLevel; // limiter level (amplitude) to reduce time and frequency localised interference

  public:

   int WrapMask;   // wrap mask for buffer addressing

   Type *InpTap;      // input buffer for analysis window
   int InpTapPtr;

   Type *OutTap;      // output buffer for reconstruction window
   int OutTapPtr;

   Type *WindowShape; // analysis/reconstruction window shape

   int SliceSepar;  // time separation between analysis/reconstruction slices

   r2FFT< Cmpx<Type> > FFT; // FFT engine
   Cmpx<Type> *FFT_Buff;    // FFT buffer

   int  SpectraLen;      // number of spectral points after FFT
   Cmpx<Type> *Spectra[2];  // spectra buffer

   Type *Output;            // (final) output buffer after pulse limiter

   Type *Energy;            // energy vs frequency

   BoxFilter<Type> Filter;  // spectra energy averaging filter

  public:
   MFSK_InputProcessor()
   { Init();
     Default(); }

   ~MFSK_InputProcessor()
   { Free(); }
   
   void Init(void)
     { InpTap=0;
       OutTap=0;
	   WindowShape=0;
	   FFT_Buff=0;
       Spectra[0]=0;
	   Spectra[1]=0;
	   Output=0;
       Energy=0; }

   void Free(void)
     { free(InpTap); InpTap=0;
       free(OutTap); OutTap=0;
	   free(WindowShape); WindowShape=0;
	   free(FFT_Buff); FFT_Buff=0;
	   free(Spectra[0]); Spectra[0]=0;
       free(Spectra[1]); Spectra[1]=0;
	   free(Output); Output=0;
	   free(Energy); Energy=0;
       FFT.Free();
       Filter.Free(); }

   void Default(void)
     { WindowLen=8192;
	   LimiterLevel=2.5; }
     
   int Preset(void)
     { int Idx;

       WrapMask=WindowLen-1;

       Type ShapeScale=2.0/WindowLen;

       if(ReallocArray(&InpTap,WindowLen)<0) goto Error;
       ClearArray(InpTap,WindowLen);
       InpTapPtr=0;
       if(ReallocArray(&OutTap,WindowLen)<0) goto Error;
       ClearArray(OutTap,WindowLen);
       OutTapPtr=0;

       if(FFT.Preset(WindowLen)<0) goto Error;
       if(ReallocArray(&FFT_Buff,WindowLen)<0) goto Error;
       SliceSepar=WindowLen/2;

       if(ReallocArray(&WindowShape,WindowLen)<0) goto Error;
       for(Idx=0; Idx<WindowLen; Idx++)
         WindowShape[Idx]=ShapeScale*sqrt(1.0-FFT.Twiddle[Idx].Re);

       SpectraLen=WindowLen/2;
       if(ReallocArray(&Spectra[0],SpectraLen)<0) goto Error;
       if(ReallocArray(&Spectra[1],SpectraLen)<0) goto Error;

       if(ReallocArray(&Output,WindowLen)<0) goto Error;
       ClearArray(Output,WindowLen);

       if(ReallocArray(&Energy,SpectraLen)<0) goto Error;

       Filter.Len=WindowLen/16;
       if(Filter.Preset()<0) goto Error;

       return 0;
       
       Error: Free(); return -1; }

   void Reset(void)
     { ClearArray(InpTap,WindowLen);
       InpTapPtr=0;
       ClearArray(OutTap,WindowLen);
       OutTapPtr=0; }

   void LimitSpectraPeaks(Cmpx<Type> *Spectra, int BoxLen=64)
     { Filter.Len=BoxLen;
       Filter.Preset();

       int MaxFreq = 3*(SpectraLen/4);
       int Freq,Idx;

       for(Freq=0; Freq<BoxLen; Freq++)
         Filter.Process(Energy[Freq]);

       Type Threshold=LimiterLevel*LimiterLevel;
       for(Idx=BoxLen/2; Freq<MaxFreq; Freq++,Idx++)
       { Filter.Process(Energy[Freq]);
         Type Signal = Energy[Idx];
         Type Limit=(Filter.Output/BoxLen)*Threshold;
         if(Signal>Limit)
         { Spectra[Idx]*=sqrt(Limit/Signal);
           Energy[Idx]=Limit; }
       }

     }

   void LimitOutputPeaks(void)
     {
       int Idx;
       Type RMS=0;
       for(Idx=0; Idx<WindowLen; Idx++)
       { Type Signal=Output[Idx];
         RMS+=Signal*Signal; }
       RMS=sqrt(RMS/WindowLen);
       Type Limit=RMS*LimiterLevel;

       for(Idx=0; Idx<WindowLen; Idx++)
       { Type Signal=Output[Idx];
         if(Signal>Limit)
           Output[Idx]=Limit;
         else if(Signal<(-Limit))
           Output[Idx]=(-Limit);
       }

     }

   void AverageEnergy(int Len=32)
     { Filter.Len=Len;
       Filter.Preset();

       int MaxFreq = 3*(SpectraLen/4);
       Type Scale=1.0/Len;
       int Len2=Len/2;
       int Idx,Freq;

       for(Freq=0; Freq<Len; Freq++)
         Filter.Process(Energy[Freq]);

       for(Idx=0; Idx<Len2; Idx++)
         Energy[Idx]=Filter.Output*Scale;

       for(      ; Freq<MaxFreq; Freq++,Idx++)
       { Filter.Process(Energy[Freq]);
         Energy[Idx]=Filter.Output*Scale; }

       for(      ; Idx<SpectraLen; Idx++)
         Energy[Idx]=Filter.Output*Scale;

     }

   // here we process the spectral data
   void ProcessSpectra(Cmpx<Type> *Spectra)
     { int Freq;

       for(Freq=0; Freq<SpectraLen; Freq++)
       { Energy[Freq]=Spectra[Freq].Energy(); }

       LimitSpectraPeaks(Spectra,WindowLen/64);
       LimitSpectraPeaks(Spectra,WindowLen/64);
       LimitSpectraPeaks(Spectra,WindowLen/64);

       AverageEnergy(WindowLen/96);
       AverageEnergy(WindowLen/64);

       for(Freq=0; Freq<SpectraLen; Freq++)
       { Type Corr=Energy[Freq];
         if(Corr<=0) continue;
         Corr=1.0/sqrt(Corr);
         Spectra[Freq]*=Corr; }

     }

   template <class InpType>
    void ProcessInpTap(InpType *Input)
     { int InpIdx;
       for(InpIdx=0; InpIdx<SliceSepar; InpIdx++)
       { InpTap[InpTapPtr]=Input[InpIdx];
         InpTapPtr+=1; InpTapPtr&=WrapMask; }
     }

   void ProcessInpTap()
     { int InpIdx;
       for(InpIdx=0; InpIdx<SliceSepar; InpIdx++)
       { InpTap[InpTapPtr]=0;
         InpTapPtr+=1; InpTapPtr&=WrapMask; }
     }

   void ProcessInpWindow_Re(void)
     { int Time;
       for(Time=0; Time<WindowLen; Time++)
       { FFT_Buff[Time].Re=InpTap[InpTapPtr]*WindowShape[Time];
         InpTapPtr+=1; InpTapPtr&=WrapMask; }
     }

   void ProcessInpWindow_Im(void)
     { int Time;
       for(Time=0; Time<WindowLen; Time++)
       { FFT_Buff[Time].Im=InpTap[InpTapPtr]*WindowShape[Time];
         InpTapPtr+=1; InpTapPtr&=WrapMask; }
     }

   void ProcessOutWindow_Re(void)
     { int Time;
       for(Time=0; Time<WindowLen; Time++)
       { OutTap[OutTapPtr]+=FFT_Buff[Time].Re*WindowShape[Time];
         OutTapPtr+=1; OutTapPtr&=WrapMask; }
     }

   void ProcessOutWindow_Im(void)
     { int Time;
       for(Time=0; Time<WindowLen; Time++)
       { OutTap[OutTapPtr]+=FFT_Buff[Time].Im*WindowShape[Time];
         OutTapPtr+=1; OutTapPtr&=WrapMask; }
     }

   void ProcessOutTap(Type *Output)
     { int OutIdx;
       for(OutIdx=0; OutIdx<SliceSepar; OutIdx++)
       { Output[OutIdx]=OutTap[OutTapPtr];
         OutTap[OutTapPtr]=0;
         OutTapPtr+=1; OutTapPtr&=WrapMask; }
     }

   template <class InpType>
    int Process(InpType *Input)
     {
       if(Input) ProcessInpTap(Input);
            else ProcessInpTap();
       ProcessInpWindow_Re();
       if(Input) ProcessInpTap(Input+SliceSepar);
            else ProcessInpTap();
       ProcessInpWindow_Im();

       FFT.Process(FFT_Buff);
       FFT.SeparTwoReals(FFT_Buff, Spectra[0], Spectra[1]);

       ProcessSpectra(Spectra[0]);
       ProcessSpectra(Spectra[1]);

       FFT.JoinTwoReals(Spectra[0], Spectra[1], FFT_Buff);
       FFT.Process(FFT_Buff);

       ProcessOutWindow_Re();
       ProcessOutTap(Output);
       ProcessOutWindow_Im();
       ProcessOutTap(Output+SliceSepar);

       LimitOutputPeaks();
       LimitOutputPeaks();

       return WindowLen; }

   // get output as 16-bit signed data
   int GetOutput(int16_t *Buffer)
     { const Type Scale=32768.0;
       const int32_t Limit=0x7FFF;
       int Idx;

       for(Idx=0; Idx<WindowLen; Idx++)
       { Type Ampl=Output[Idx];
         Ampl*=Scale;
         int32_t Out=(int32_t)floor(Ampl+0.5);
         if(Out>Limit) Out=Limit;
         else if(Out<(-Limit)) Out=(-Limit);
         Buffer[Idx]=(int16_t)Out; }

       return WindowLen; }

} ;

// =====================================================================

// front-end spectral analysis and soft decoder for MFSK
template <class Type=float>
 class MFSK_Demodulator
{ public:

   MFSK_Parameters<Type> *Parameters;

  public:

   int InputLen;       // input must be provided in batches of that length [samples]

  private:

   int SymbolSepar;
   int SymbolLen;
   int SpectraPerSymbol;

   int DecodeMargin;   // frequency margin for decoding the signal [FFT bins]
   int DecodeWidth;    // Spectra width                            [FFT bins]

   int SliceSepar;     // time separation between samples          [samples]

   int WrapMask;                     // wrapping mask for the FFT window

   Type *InpTap;                        // input buffer
   int InpTapPtr;

   Type *SymbolShape;                   // the shape of the symbol and the FFT window

   r2FFT< Cmpx<Type> > FFT;             // FFT engine
   Cmpx<Type> *FFT_Buff;                // FFT buffer

   int  SpectraLen;                   // number of spectral points per FFT
   Cmpx<Type> *Spectra[2];              // 2 buffers for FFT spectra

   CircularBuffer<Type, int> History;         // Spectra history

  public:

   MFSK_Demodulator()
   { Init(); }

   ~MFSK_Demodulator()
   { Free(); }

   void Init(void)
     { InpTap=0;
       SymbolShape=0;
       FFT_Buff=0;
       Spectra[0]=0;
       Spectra[1]=0; }

   void Free(void)
     { free(InpTap); InpTap=0;
       free(SymbolShape); SymbolShape=0;
       free(FFT_Buff); FFT_Buff=0;
       free(Spectra[0]); Spectra[0]=0;
       free(Spectra[1]); Spectra[1]=0;
       FFT.Free();
       History.Free(); }

   int Preset(MFSK_Parameters<Type> *NewParameters)
     {
       Parameters=NewParameters;

       SymbolSepar=Parameters->SymbolSepar;
       SymbolLen=Parameters->SymbolLen;
       SpectraPerSymbol=Parameters->SpectraPerSymbol;

       InputLen=SymbolSepar;
       DecodeMargin=Parameters->RxSyncMargin*Parameters->CarrierSepar;

       WrapMask=SymbolLen-1;

       Type ShapeScale=1.0/SymbolLen;

       if(ReallocArray(&InpTap,SymbolLen)<0) goto Error;
       ClearArray(InpTap,SymbolLen);
       InpTapPtr=0;

       if(FFT.Preset(SymbolLen)<0) goto Error;
       if(ReallocArray(&FFT_Buff,SymbolLen)<0) goto Error;
       SliceSepar=SymbolSepar/SpectraPerSymbol;

       if(ReallocArray(&SymbolShape,SymbolLen)<0) goto Error;

       { int Time;
         double Ampl=MFSK_SymbolFreqShape[0];
         for(Time=0; Time<SymbolLen; Time++)
           SymbolShape[Time]=Ampl;
       }
       int Freq;
	   for(Freq=1; Freq<MFSK_SymbolFreqShapeLen; Freq++)
       { int Time;
         double Ampl=MFSK_SymbolFreqShape[Freq];
         if(Freq&1) Ampl=(-Ampl);
         int Phase=0;
         for(Time=0; Time<SymbolLen; Time++)
         { SymbolShape[Time]+=Ampl*FFT.Twiddle[Phase].Re;
           Phase+=Freq; if(Phase>=SymbolLen) Phase-=SymbolLen; }
       }
       { int Time;
         for(Time=0; Time<SymbolLen; Time++)
           SymbolShape[Time]*=ShapeScale;
       }

       SpectraLen=SymbolLen/2;
       if(ReallocArray(&Spectra[0],SpectraLen)<0) goto Error;
       if(ReallocArray(&Spectra[1],SpectraLen)<0) goto Error;

       DecodeWidth=((Parameters->Carriers-1)*Parameters->CarrierSepar+1) + 2*DecodeMargin;

       History.Len=(Parameters->RxSyncIntegLen+2)*Parameters->SpectraPerBlock;
       History.Width=DecodeWidth;
       if(History.Preset()<0) goto Error;
       History.Clear();

       return 0;

       Error: Free(); return -1; }

   void Reset(void)
     { History.Clear(); }

   Type *HistoryPtr(int Idx)
     { return History.OffsetPtr(Idx); }

   template <class InpType>
    int SlideOneSlice(InpType *Input)
     { int InpIdx;
       for(InpIdx=0; InpIdx<SliceSepar; InpIdx++)
       { InpTap[InpTapPtr]=Input[InpIdx];
         InpTapPtr+=1; InpTapPtr&=WrapMask; }
	   return SliceSepar; }

   template <class InpType>
    void Process(InpType *Input)
     { int InpIdx,Time,Slice;

       for(InpIdx=0, Slice=0; Slice<SpectraPerSymbol; Slice+=2)
       {
         InpIdx+=SlideOneSlice(Input+InpIdx);

         for(Time=0; Time<SymbolLen; Time++)
         { FFT_Buff[Time].Re=InpTap[InpTapPtr]*SymbolShape[Time];
           InpTapPtr+=1; InpTapPtr&=WrapMask; }

         InpIdx+=SlideOneSlice(Input+InpIdx);

         for(Time=0; Time<SymbolLen; Time++)
         { FFT_Buff[Time].Im=InpTap[InpTapPtr]*SymbolShape[Time];
           InpTapPtr+=1; InpTapPtr&=WrapMask; }

         FFT.Process(FFT_Buff);
         FFT.SeparTwoReals(FFT_Buff, Spectra[0], Spectra[1]);

         Type *Data0 = History.OffsetPtr(0);
         Type *Data1 = History.OffsetPtr(1);

         int Idx;
         int Freq=Parameters->FirstCarrier-DecodeMargin;
         for(Idx=0; Idx<DecodeWidth; Idx++, Freq++)
         { Data0[Idx]=Spectra[0][Freq].Energy();
           Data1[Idx]=Spectra[1][Freq].Energy(); }

         History+=2;
       }

     }

   template <class OutType>
    int PickBlock(OutType *Spectra, int TimeOffset, int FreqOffset)
     { int SpectraPerBlock=Parameters->SpectraPerBlock;
 	   if((TimeOffset>(-SpectraPerBlock))||((-TimeOffset)>(int)History.Len)) return -1;

       int Carriers=Parameters->Carriers;
       int CarrierSepar=Parameters->CarrierSepar;

       if((FreqOffset<0)||((FreqOffset+(Carriers-1)*CarrierSepar)>=DecodeWidth)) return -1;

       int SymbolsPerBlock=Parameters->SymbolsPerBlock;
	   int Symbol;
       for(Symbol=0; Symbol<SymbolsPerBlock; Symbol++, TimeOffset+=SpectraPerSymbol)
       { Type *Hist=History.OffsetPtr(TimeOffset)+FreqOffset;
	     int Freq;
         for(Freq=0; Freq<Carriers; Freq++,Hist+=CarrierSepar)
		   (*Spectra++)=(*Hist);
	   }

	   return 0; }

} ;

// =====================================================================

template <class Type>
 void PrintBinary(Type Number, int Bits)
{ Type Mask=1; Mask<<=(Bits-1);
  for( ; Bits; Bits--)
  { printf("%c",Number&Mask ? '1':'0');
    Mask>>=1; }
}

// =====================================================================

// FEC block encoder
class MFSK_Encoder
{ public:

   // parameters to be set before calling Preset()
   int BitsPerSymbol;    // number of bits per MFSK symbol (default is 5 bits thus 32 possible symbols)
   int BitsPerCharacter; // number of bits per transmitted character (default is 7 bits for ASCII code)

  public:

   int Symbols;          // number of possible MFSK symbols
   int SymbolsPerBlock;  // number of MFSK symbols per FEC block

  private:

   static const uint64_t ScramblingCode = 0xE257E6D0291574ECLL;

   int8_t *FHT_Buffer;    // temporary buffer for (inverse) Fast Hadamard Transform

  public:

   uint8_t *OutputBlock;  // encoded block is stored here

  public:

   MFSK_Encoder()
     { Init();
       Default(); }

   ~MFSK_Encoder()
     { Free(); }

   void Default(void)
     { BitsPerSymbol=5;
       BitsPerCharacter=7; }

   void Init(void)
     { FHT_Buffer=0;
       OutputBlock=0; }

   void Free(void)
     { free(FHT_Buffer); FHT_Buffer=0;
       free(OutputBlock); OutputBlock=0; }

   int Preset(void)
     { Symbols = 1<<BitsPerSymbol;
       SymbolsPerBlock = Exp2(BitsPerCharacter-1);
       if(ReallocArray(&FHT_Buffer,SymbolsPerBlock)<0) goto Error;
       if(ReallocArray(&OutputBlock,SymbolsPerBlock)<0) goto Error;
       return 0;
       Error: Free(); return -1; }

   // encode a single character with the IFHT
   void EncodeCharacter(uint8_t Char)
     { int TimeBit;

       uint8_t Mask=(SymbolsPerBlock<<1)-1;
       Char&=Mask;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
         FHT_Buffer[TimeBit]=0;
       if(Char<SymbolsPerBlock) FHT_Buffer[Char]=1;
                   else FHT_Buffer[Char-SymbolsPerBlock]=(-1);
       IFHT(FHT_Buffer, SymbolsPerBlock);
     }

   // scramble the codeword (of a single character) with the scrambling code
   void ScrambleFHT(int CodeOffset=0)
     { int TimeBit;
       int CodeWrap=(SymbolsPerBlock-1);
       int CodeBit=CodeOffset&CodeWrap;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { uint64_t CodeMask=1; CodeMask<<=CodeBit;
         if(ScramblingCode&CodeMask)
           FHT_Buffer[TimeBit] = (-FHT_Buffer[TimeBit]);
       CodeBit+=1; CodeBit&=CodeWrap; }
     }

   // encode a block of SymbolsPerBlock characters
   void EncodeBlock(uint8_t *InputBlock)
     { int FreqBit,TimeBit;

       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { OutputBlock[TimeBit]=0; }

       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
       { EncodeCharacter(InputBlock[FreqBit]);
         ScrambleFHT(FreqBit*13);
         int Rotate=0;
		 for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
         { if(FHT_Buffer[TimeBit]<0)
           { int Bit=FreqBit+Rotate; if(Bit>=BitsPerSymbol) Bit-=BitsPerSymbol;
             uint8_t Mask=1; Mask<<=Bit;
             OutputBlock[TimeBit]|=Mask; }
           Rotate+=1; if(Rotate>=BitsPerSymbol) Rotate-=BitsPerSymbol; }
       }

     }

   // print the encoded block (for debug only)
   void PrintOutputBlock(void)
     { int TimeBit;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { printf("%2d: ",TimeBit);
	     PrintBinary(OutputBlock[TimeBit],BitsPerSymbol);
		 printf("\n"); }
     }
   
} ;

/*
// Hard FEC decoder
class MFSK_HardDecoder
{ public:

   // parameters to be set before calling Preset()
   int BitsPerSymbol;    // number of bits per MFSK symbol (default is 5 bits thus 32 possible symbols)
   int BitsPerCharacter; // number of bits per transmitted character (default is 7 bits for ASCII code)

  public:

   int Symbols;          // number of possible MFSK symbols
   int SymbolsPerBlock;  // number of MFSK symbols per FEC block

  private:
  
   static const uint64_t ScramblingCode = 0xE257E6D0291574ECLL;

   uint8_t *InputBuffer;
   int InputPtr;
   int InputWrap;

   int8_t *FHT_Buffer;

  public:

   float Signal,NoiseEnergy;
   uint8_t *OutputBlock;

  public:

   MFSK_HardDecoder()
     { Init();
       Default(); }

   ~MFSK_HardDecoder()
     { Free(); }

   void Default(void)
     { BitsPerSymbol=5;
       BitsPerCharacter=7; }

   void Init(void)
     { InputBuffer=0;
       FHT_Buffer=0;
       OutputBlock=0; }

   void Free(void)
     { free(InputBuffer); InputBuffer=0;
       free(FHT_Buffer); FHT_Buffer=0;
       free(OutputBlock); OutputBlock=0; }

   void Reset(void)
     { int Idx;
       for(Idx=0; Idx<SymbolsPerBlock; Idx++)
         InputBuffer[Idx]=0;
       InputPtr=0; }

   int Preset(void)
     { Symbols = 1<<BitsPerSymbol;
       SymbolsPerBlock = Exp2(BitsPerCharacter-1);
       if(ReallocArray(&InputBuffer,SymbolsPerBlock)<0) goto Error;
       if(ReallocArray(&FHT_Buffer,SymbolsPerBlock)<0) goto Error;
       if(ReallocArray(&OutputBlock,BitsPerSymbol)<0) goto Error;

       InputWrap = SymbolsPerBlock-1;
       Reset();

       return 0;
       Error: Free(); return -1; }

   // input a new MFSK symbol
   void Input(uint8_t Symbol)
     { InputBuffer[InputPtr]=Symbol;
       InputPtr+=1; InputPtr&=InputWrap; }

   // decode given character in the block
   void DecodeCharacter(int FreqBit)
     { int TimeBit;

       int Ptr=InputPtr;
       int Rotate=FreqBit;
       int CodeWrap=(SymbolsPerBlock-1);
       int CodeBit=FreqBit*13; CodeBit&=CodeWrap;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { uint8_t Bit=InputBuffer[Ptr]; Bit>>=Rotate; Bit&=1;
         uint64_t CodeMask=1; CodeMask<<=CodeBit;
         if(ScramblingCode&CodeMask) Bit^=1;
         FHT_Buffer[TimeBit]= Bit ? -1:1;
         CodeBit+=1; CodeBit&=CodeWrap;
         Rotate+=1; if(Rotate>=BitsPerSymbol) Rotate-=BitsPerSymbol;
         Ptr+=1; Ptr&=InputWrap; }

       FHT(FHT_Buffer,SymbolsPerBlock);
       int32_t Peak=0;
       int PeakPos=0;
       int32_t SqrSum=0;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { int32_t Signal=FHT_Buffer[TimeBit];
         SqrSum+=Signal*Signal;
         if(abs(Signal)>abs(Peak))
         { Peak=Signal; PeakPos=TimeBit; }
       }

       uint8_t Char=PeakPos;
       if(Peak<0) Char+=SymbolsPerBlock;
       SqrSum-=Peak*Peak;

       OutputBlock[FreqBit]=Char;
       NoiseEnergy+=(float)SqrSum/(SymbolsPerBlock-1);
       Signal+=abs(Peak);
     }

   // decode all characters in a block, sum up the noise and signal
   void Process(void)
     { int FreqBit;
       Signal=0; NoiseEnergy=0;
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
         DecodeCharacter(FreqBit);
       Signal/=BitsPerSymbol;
       NoiseEnergy/=BitsPerSymbol;
     }

   // get the decoded block
   int Output(uint8_t *Buffer)
     { int FreqBit;
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
         Buffer[FreqBit]=OutputBlock[FreqBit];
       return BitsPerSymbol; }

   // get the decoded block in a packed form
   int Output(uint64_t &PackedBuffer)
     { int FreqBit;
       PackedBuffer=0;
       for(FreqBit=BitsPerSymbol; FreqBit>0; )
       { PackedBuffer<<=8;
         FreqBit--;
         PackedBuffer|=OutputBlock[FreqBit]; }
       return BitsPerSymbol; }

   // like above, but pointer instead of reference
   int Output(uint64_t *PackedBuffer)
     { return Output(*PackedBuffer); }

  // print decoded block (for debug)
   void PrintOutputBlock(FILE *File=stdout)
     { int FreqBit;
       fprintf(File,"'");
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
       { uint8_t Char=OutputBlock[FreqBit];
         fprintf(File,"%c", (Char>=' ')&&(Char<127) ? Char:' '); }
       fprintf(File,"', S/N = %5.1f/%4.1f",Signal,sqrt(NoiseEnergy));
       if(NoiseEnergy>0) fprintf(File," = %5.1f",Signal/sqrt(NoiseEnergy));
       fprintf(File,"\n"); }

   // print the input buffer that is MFSK symbols received with Input()
   void PrintInputBuffer(void)
     { int TimeBit;
       int Ptr=InputPtr;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { printf("%2d: ",TimeBit);
	     PrintBinary(InputBuffer[Ptr],BitsPerSymbol);
		 printf("\n");
		 Ptr+=1; Ptr&=InputWrap; }
     }
   
} ;
*/
// Soft FEC decoder (calls are similar to the MFSK_HardDecoder class)
template <class InpType=float, class CalcType=float>
 class MFSK_SoftDecoder
{ public:

   MFSK_Parameters<CalcType> *Parameters;

  private:

   int BitsPerSymbol;
   int SymbolsPerBlock;
   int SpectraPerSymbol;

   int InputBufferLen;
   InpType *InputBuffer;
   int InputPtr;

   CalcType *FHT_Buffer;

  public:
   CalcType Signal,NoiseEnergy;
   uint8_t *OutputBlock;

  public:

   MFSK_SoftDecoder()
     { Init(); }

   ~MFSK_SoftDecoder()
     { Free(); }

   void Init(void)
     { InputBuffer=0;
       FHT_Buffer=0;
       OutputBlock=0; }

   void Free(void)
     { free(InputBuffer); InputBuffer=0;
       free(FHT_Buffer); FHT_Buffer=0;
       free(OutputBlock); OutputBlock=0; }

   void Reset(void)
     { int Idx;
       for(Idx=0; Idx<InputBufferLen; Idx++)
         InputBuffer[Idx]=0;
       InputPtr=0; }

   int Preset(MFSK_Parameters<float> *NewParameters)
     { Parameters=NewParameters;

	   BitsPerSymbol=Parameters->BitsPerSymbol;
	   SymbolsPerBlock=Parameters->SymbolsPerBlock;
	   SpectraPerSymbol=Parameters->SpectraPerSymbol;
	   InputBufferLen=SymbolsPerBlock*SpectraPerSymbol*BitsPerSymbol;
       if(ReallocArray(&InputBuffer,InputBufferLen)<0) goto Error;
       if(ReallocArray(&FHT_Buffer,SymbolsPerBlock)<0) goto Error;
       if(ReallocArray(&OutputBlock,BitsPerSymbol)<0) goto Error;
       Reset();

       return 0;
       Error: Free(); return -1; }

   void SpectralInput(InpType *SpectraEnergy)
     { MFSK_SoftDemodulate(InputBuffer+InputPtr,SpectraEnergy, BitsPerSymbol,
	      Parameters->CarrierSepar, Parameters->UseGrayCode, Parameters->RxSyncSquareEnergy);
       InputPtr+=BitsPerSymbol;
       if(InputPtr>=InputBufferLen) InputPtr-=InputBufferLen; }

   void Input(InpType *Symbol)
     { int FreqBit;
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
       { InputBuffer[InputPtr]=Symbol[FreqBit];
         InputPtr+=1; }
       if(InputPtr>=InputBufferLen) InputPtr-=InputBufferLen; }

   void DecodeCharacter(int FreqBit)
     { int TimeBit;

       int Ptr=InputPtr;
       int Rotate=FreqBit;
       int CodeWrap=(SymbolsPerBlock-1);
       int CodeBit=FreqBit*13; CodeBit&=CodeWrap;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { InpType Bit=InputBuffer[Ptr+Rotate];
         uint64_t CodeMask=1; CodeMask<<=CodeBit;
         if(Parameters->ScramblingCode&CodeMask) Bit=(-Bit);
         FHT_Buffer[TimeBit]=Bit;
         CodeBit+=1; CodeBit&=CodeWrap;
         Rotate+=1; if(Rotate>=BitsPerSymbol) Rotate-=BitsPerSymbol;
         Ptr+=(BitsPerSymbol*SpectraPerSymbol);
         if(Ptr>=InputBufferLen) Ptr-=InputBufferLen; }

       FHT(FHT_Buffer,SymbolsPerBlock);
       CalcType Peak=0;
       int PeakPos=0;
       CalcType SqrSum=0;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { CalcType Signal=FHT_Buffer[TimeBit];
         SqrSum+=Signal*Signal;
         if(fabs(Signal)>fabs(Peak))
         { Peak=Signal; PeakPos=TimeBit; }
       }

       uint8_t Char=PeakPos;
       if(Peak<0) Char+=SymbolsPerBlock;
       SqrSum-=Peak*Peak;

       OutputBlock[FreqBit]=Char;
       NoiseEnergy+=(float)SqrSum/(SymbolsPerBlock-1);
       Signal+=fabs(Peak);
     }

   void Process(void)
     { int FreqBit;
       Signal=0; NoiseEnergy=0;
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
         DecodeCharacter(FreqBit);
       Signal/=BitsPerSymbol;
       NoiseEnergy/=BitsPerSymbol;
     }

   int Output(uint8_t *Buffer)
     { int FreqBit;
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
         Buffer[FreqBit]=OutputBlock[FreqBit];
       return BitsPerSymbol; }

   void PrintOutputBlock(FILE *File=stdout)
     { int FreqBit;
       fprintf(File,"'");
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
       { uint8_t Char=OutputBlock[FreqBit];
         fprintf(File,"%c", (Char>=' ')&&(Char<127) ? Char:' '); }
       fprintf(File,"', S/N = %5.1f/%4.1f",Signal,sqrt(NoiseEnergy));
       if(NoiseEnergy>0) fprintf(File," = %5.1f",Signal/sqrt(NoiseEnergy));
       fprintf(File,"\n"); }

} ;

// Soft but iterative (!) FEC decoder
template <class Type>
 class MFSK_SoftIterDecoder
{ public:

   MFSK_Parameters<Type> *Parameters;

   Type *Input;				// demodulated spectra energies / tone probabilities

  private:

   int BitsPerSymbol;
   int BitsPerCharacter;
   int Symbols;
   int SymbolsPerBlock;

   Type *InputExtrinsic;	// extrinsic input information fed back from the decoder
   Type *FHT_Codeword;		// FHT codewords to be decoded by FHT

  public:

   Type Input_SignalEnergy;
   Type Input_NoiseEnergy;
   Type FEC_SignalEnergy;
   Type FEC_NoiseEnergy;

   uint8_t *OutputBlock;

  public:

   MFSK_SoftIterDecoder()
     { Init(); }

   ~MFSK_SoftIterDecoder()
     { Free(); }

   void Init(void)
     { Input=0;
	   InputExtrinsic=0;
       FHT_Codeword=0;
       OutputBlock=0; }

   void Free(void)
     { free(Input); Input=0;
       free(InputExtrinsic); InputExtrinsic=0;
       free(FHT_Codeword); FHT_Codeword=0;
       free(OutputBlock); OutputBlock=0; }

   int Preset(MFSK_Parameters<Type> *NewParameters)
     { Parameters=NewParameters;
       BitsPerSymbol=Parameters->BitsPerSymbol;
       BitsPerCharacter=Parameters->BitsPerCharacter;
	   Symbols = Parameters->Carriers;
       SymbolsPerBlock = Parameters->SymbolsPerBlock;
       if(ReallocArray(&Input,SymbolsPerBlock*Symbols)<0) goto Error;
       if(ReallocArray(&InputExtrinsic,SymbolsPerBlock*Symbols)<0) goto Error;
       if(ReallocArray(&FHT_Codeword,SymbolsPerBlock*BitsPerSymbol)<0) goto Error;
       if(ReallocArray(&OutputBlock,BitsPerSymbol)<0) goto Error;

       return 0;
       Error: Free(); return -1; }

   void SimulateInput(uint8_t *InputBlock, Type SNR=1.0, Type DeadCarrierSNR=0.0)
     { Type NoiseRMS=1.0;
	   Type Signal=SNR*NoiseRMS*sqrt(2.0*Symbols);
       Type DeadCarrier=DeadCarrierSNR*NoiseRMS*sqrt(2.0*Symbols);
	   int DeadCarrierFreq=rand()&(Symbols-1);

	   int Symbol;
       int Idx;
       for(Symbol=0,Idx=0; Symbol<SymbolsPerBlock; Symbol++,Idx+=Symbols)
       { int Freq;
	     int SymbolFreq=InputBlock[Symbol];
         if(Parameters->UseGrayCode) SymbolFreq=GrayCode(SymbolFreq);
         for(Freq=0; Freq<Symbols; Freq++)
         { Cmpx<Type> Noise;
		   WhiteNoise(Noise,NoiseRMS);
		   if(Freq==SymbolFreq) Noise.Re+=Signal;
           if(Freq==DeadCarrierFreq) Noise.Im+=DeadCarrier;
           Type Energy=Noise.Mag2();
		   Input[Idx+Freq]=Energy*Energy;
		 }
       }

	 }

   int NormalizeSum(Type *Data, int Len, Type Norm=1.0)
     { int Idx;
	   Type Sum=0;
	   for(Idx=0; Idx<Len; Idx++)
	     Sum+=Data[Idx];
	   if(Sum<=0) return -1;
	   Type Corr=Norm/Sum;
	   for(Idx=0; Idx<Len; Idx++)
	     Data[Idx]*=Corr;
	   return 0; }

   int NormalizeAbsSum(Type *Data, int Len, Type Norm=1.0)
     { int Idx;
	   Type Sum=0;
	   for(Idx=0; Idx<Len; Idx++)
	     Sum+=fabs(Data[Idx]);
	   if(Sum<=0) return -1;
	   Type Corr=Norm/Sum;
	   for(Idx=0; Idx<Len; Idx++)
	     Data[Idx]*=Corr;
	   return 0; }

   int ThirdPower(Type *Data, int Len)
     { int Idx;
	   for(Idx=0; Idx<Len; Idx++)
	   { Type Square=Data[Idx];
	     // Square=Square*Square;
		 Square=fabs(Square);		// <= works better (in simulations)
		 Data[Idx]*=Square; }
	   return 0; }
/*
   void CalculateSymbolBitProb(Type *SymbolBit, Type *FreqProb)
     { int Bit;
	   for(Bit=0; Bit<BitsPerSymbol; Bit++)
	     SymbolBit[Bit]=0;
	   int Freq;
	   for(Freq=0; Freq<Symbols; Freq++)
	   { Type Prob=FreqProb[Freq];
	     uint8_t Symbol=Freq;
	     if(Parameters->UseGrayCode) Symbol=BinaryCode(Symbol);
         uint8_t Mask=1;
         for(Bit=0; Bit<BitsPerSymbol; Bit++)
		 { if(Symbol&Mask) SymbolBit[Bit]-=Prob;
		              else SymbolBit[Bit]+=Prob;
		   Mask<<=1; }
	   }
	 }

   void CalculateFreqProb( Type *FreqProb, Type *SymbolBit)
     { int Freq;
	   for(Freq=0; Freq<Symbols; Freq++)
       { uint8_t Symbol=Freq;
	     if(Parameters->UseGrayCode) Symbol=BinaryCode(Symbol);
		 Type Prob=1;
		 int Bit;
         uint8_t Mask=1;
         for(Bit=0; Bit<BitsPerSymbol; Bit++)
		 { Type BitProb=1.0;
		   if(Symbol&Mask) BitProb-=SymbolBit[Bit];
		              else BitProb+=SymbolBit[Bit];
           Prob*=(BitProb/2);
		   Mask<<=1; }
	     FreqProb[Freq]=Prob;
	   }
	 }
*/
   void ScrambleCodeword(Type *CodeWord, int ScrambleIdx)
     { int Idx;
	   int CodeWrap=(SymbolsPerBlock-1);
       ScrambleIdx&=CodeWrap;
       for(Idx=0; Idx<SymbolsPerBlock; Idx++)
	   { uint64_t CodeMask=1; CodeMask<<=ScrambleIdx;
	     if(Parameters->ScramblingCode&CodeMask) CodeWord[Idx]=(-CodeWord[Idx]);
		 ScrambleIdx+=1; ScrambleIdx&=CodeWrap; }
	 }

   uint8_t DecodeChar(Type *FHT_Buffer)
     { int TimeBit;
	   Type Peak=0;
       int PeakPos=0;
       Type NoiseEnergy=0;
       for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { Type Signal=FHT_Buffer[TimeBit];
         NoiseEnergy+=Signal*Signal;
         if(fabs(Signal)>fabs(Peak))
         { Peak=Signal; PeakPos=TimeBit; }
       }
       uint8_t Char=PeakPos;
       if(Peak<0) Char+=SymbolsPerBlock;
       Type SignalEnergy=Peak*Peak;
       NoiseEnergy-=SignalEnergy;
       SignalEnergy-=NoiseEnergy/(SymbolsPerBlock-1);
       NoiseEnergy*=(Type)SymbolsPerBlock/(SymbolsPerBlock-1);
/*
       printf("  %+5.1f/%3.1f=>%02X/%c",
	          Peak, sqrt(NoiseEnergy/SymbolsPerBlock), Char, Char>' ' ? Char:' ');
*/

       FEC_SignalEnergy+=SignalEnergy;
       FEC_NoiseEnergy+=NoiseEnergy;

	   return Char; }

   template <class DstType, class SrcType>
    void Copy(DstType *Dst, SrcType *Src, int Len)
	 { int Idx;
	   for(Idx=0; Idx<Len; Idx++)
	     Dst[Idx]=Src[Idx];
	 }

   template <class DstType, class SrcType>
    void Multiply(DstType *Dst, SrcType *Src, int Len)
	 { int Idx;
	   for(Idx=0; Idx<Len; Idx++)
	     Dst[Idx]*=Src[Idx];
	 }

   void Process(int MaxIter=4)
     { int TimeBit;
       int Bit;
       int Freq;
       int InpIdx;
       int BlockIdx;
       int InputSize = Symbols*SymbolsPerBlock;
       int BlockSize = BitsPerSymbol*SymbolsPerBlock;
/*
     for(TimeBit=0,InpIdx=0; TimeBit<SymbolsPerBlock; TimeBit++,InpIdx+=Symbols)
     { Copy(InputExtrinsic+InpIdx, Input+InpIdx, Symbols);
	   NormalizeSum(InputExtrinsic+InpIdx, Symbols, 1.0); }
*/
     for(InpIdx=0; InpIdx<InputSize; InpIdx++)
	   InputExtrinsic[InpIdx]=1.0/Symbols;

     for( ; MaxIter; MaxIter--)
     {

       int SquareEnergy=Parameters->DecodeSquareEnergy;
       for(InpIdx=0; InpIdx<InputSize; InpIdx++)
	   { Type InputEnergy=Input[InpIdx];
         if(SquareEnergy) InputEnergy*=InputEnergy;
	     InputExtrinsic[InpIdx]*=InputEnergy;
	   }

	   Type SymbolBit[BitsPerSymbol];
       int Rotate=0;
       for(TimeBit=0,InpIdx=0; TimeBit<SymbolsPerBlock; TimeBit++,InpIdx+=Symbols)
       { 

         MFSK_SoftDemodulate(SymbolBit, InputExtrinsic+InpIdx, BitsPerSymbol, 1, Parameters->UseGrayCode, 0);

	     // NormalizeSum(InputExtrinsic+InpIdx,Symbols,1.0);
/*
         printf("%2d:",TimeBit);
		 for(Freq=0; Freq<Symbols; Freq++)
		   printf(" %+5.2f",InputExtrinsic[InpIdx+Freq]);
		 printf("\n");
*/
	   	 // CalculateSymbolBitProb(SymbolBit,InputExtrinsic+InpIdx);
/*
         printf("%2d:",TimeBit);
		 for(Bit=0; Bit<BitsPerSymbol; Bit++)
		   printf(" %+5.2f",SymbolBit[Bit]);
		 printf("\n");
*/
         BlockIdx=TimeBit+Rotate*SymbolsPerBlock;
         for(Bit=0; Bit<BitsPerSymbol; Bit++)
		 { FHT_Codeword[BlockIdx]=SymbolBit[Bit];
		   BlockIdx+=SymbolsPerBlock; if(BlockIdx>=BlockSize) BlockIdx-=BlockSize; }

	     if(Rotate>0) Rotate-=1; else Rotate+=(BitsPerSymbol-1);
	   }

       FEC_SignalEnergy=0;
       FEC_NoiseEnergy=0;
       for(Bit=0,BlockIdx=0; Bit<BitsPerSymbol; Bit++,BlockIdx+=SymbolsPerBlock)
       { ScrambleCodeword(FHT_Codeword+BlockIdx,13*Bit);
         FHT(FHT_Codeword+BlockIdx,SymbolsPerBlock);

         uint8_t Char=DecodeChar(FHT_Codeword+BlockIdx);
		 OutputBlock[Bit]=Char;
         // printf(" %02X [%c]",Char,Char>' ' ? Char:' ');

         ThirdPower(FHT_Codeword+BlockIdx,SymbolsPerBlock);
		 NormalizeAbsSum(FHT_Codeword+BlockIdx,SymbolsPerBlock,1.0);
/*
		 printf("%d:",Bit);
         for(TimeBit=0; TimeBit<SymbolsPerBlock; TimeBit++)
		   printf(" %+5.1f",FHT_Codeword[BlockIdx+TimeBit]);
		 printf("\n");
*/
         IFHT(FHT_Codeword+BlockIdx,SymbolsPerBlock);
         ScrambleCodeword(FHT_Codeword+BlockIdx,13*Bit);
	   }

       Rotate=0;
       for(TimeBit=0,InpIdx=0; TimeBit<SymbolsPerBlock; TimeBit++,InpIdx+=Symbols)
       { BlockIdx=TimeBit+Rotate*SymbolsPerBlock;
         for(Bit=0; Bit<BitsPerSymbol; Bit++)
		 { SymbolBit[Bit]=FHT_Codeword[BlockIdx];
		   BlockIdx+=SymbolsPerBlock; if(BlockIdx>=BlockSize) BlockIdx-=BlockSize; }
/*
         printf("%2d:",TimeBit);
		 for(Bit=0; Bit<BitsPerSymbol; Bit++)
		   printf(" %+5.2f",SymbolBit[Bit]);
		 printf("\n");

         printf("%2d:",TimeBit);
		 for(Freq=0; Freq<Symbols; Freq++)
		   printf(" %+5.2f",InputExtrinsic[InpIdx+Freq]);
		 printf("\n");
*/
         // CalculateFreqProb(InputExtrinsic+InpIdx,SymbolBit);
		 MFSK_SoftModulate(InputExtrinsic+InpIdx, SymbolBit, BitsPerSymbol, Parameters->UseGrayCode);
/*
         printf("%2d:",TimeBit);
		 for(Freq=0; Freq<Symbols; Freq++)
		   printf(" %+5.2f",InputExtrinsic[InpIdx+Freq]);
		 printf("\n");
*/
	     if(Rotate>0) Rotate-=1; else Rotate+=(BitsPerSymbol-1);
	   }

       Input_SignalEnergy=0;
       Input_NoiseEnergy=0;
       for(TimeBit=0,InpIdx=0; TimeBit<SymbolsPerBlock; TimeBit++)
       { for(Freq=0; Freq<Symbols; Freq++,InpIdx++)
	     { Type Energy=Input[InpIdx];
		   Type SigProb=InputExtrinsic[InpIdx];
           Input_SignalEnergy+=SigProb*Energy;
		   Input_NoiseEnergy+=(1-SigProb)*Energy;
		 }
	   }
       Input_SignalEnergy-=Input_NoiseEnergy/(Symbols-1);
       Input_NoiseEnergy*=(Type)Symbols/(Symbols-1);

     }


	 }

   Type InputSNRdB(void)
     { return 10*log(Input_SignalEnergy/Input_NoiseEnergy)/log(10); }

   void PrintSNR(void)
     { int Bit;
       // printf("FEC: %5.1f/%4.1f", FEC_SignalEnergy, FEC_NoiseEnergy);
       // printf("FEC: %+5.1f dB", 10*log(FEC_SignalEnergy/FEC_NoiseEnergy)/log(10));
       printf("Input: %+5.1f dB", 10*log(Input_SignalEnergy/Input_NoiseEnergy)/log(10));
       printf(" : ");
       for(Bit=0; Bit<BitsPerSymbol; Bit++)
       { char Char=OutputBlock[Bit];
	     printf("%c", Char>' ' ? Char:' '); }
       printf("\n");
     }

   int WriteOutputBlock(FIFO<uint8_t, int> &Output)
     { int Bit; int Written;
       for(Written=0, Bit=0; Bit<BitsPerSymbol; Bit++)
       { uint8_t Char=OutputBlock[Bit];
         int Error=Output.Write(Char);
         if(Error<0) break;
         Written+=Error; }
       return Written; }

/*
   void PrintOutputBlock(FILE *File=stdout)
     { int FreqBit;
       fprintf(File,"'");
       for(FreqBit=0; FreqBit<BitsPerSymbol; FreqBit++)
       { uint8_t Char=OutputBlock[FreqBit];
         fprintf(File,"%c", (Char>=' ')&&(Char<127) ? Char:' '); }
       fprintf(File,"', S/N = %5.1f/%4.1f",Signal,sqrt(NoiseEnergy));
       if(NoiseEnergy>0) fprintf(File," = %5.1f",Signal/sqrt(NoiseEnergy));
       fprintf(File,"\n"); }
*/
} ;

// =====================================================================

// MFSK transmitter (FEC encoder + MFSK modulator + rate corrector)
template <class Type=float>
 class MFSK_Transmitter
{ public:

   MFSK_Parameters<Type> *Parameters;

  public:

   int MaxOutputLen;    // maximum length of the audio batch returned by Output()

  private:

   int BitsPerSymbol;
   int SymbolsPerBlock;

   static const int State_Running = 0x0001;
   static const int State_StopReq = 0x0010;
   int State;

   FIFO<uint8_t, int> Input;     // buffer(queue) for the characters to be encoded
   uint8_t InputBlock[8];        // FEC code block buffer
   FIFO<uint8_t, int> Monitor;   // buffer for monitoring the characters being sent

   MFSK_Encoder Encoder;         // FEC encoder
   int SymbolPtr;

   MFSK_Modulator<Type> Modulator; // MFSK modulator

   Type *ModulatorOutput;             // modulator output
   RateConverter<Type> RateConverter; // output rate converter
   Type *ConverterOutput;             // rate converter output

  public:

   MFSK_Transmitter()
     { Init(); }

   ~MFSK_Transmitter()
     { Free(); }

   void Init(void)
     { ModulatorOutput=0;
       ConverterOutput=0; }

   void Free(void)
     { Input.Free();
       Monitor.Free();
       Encoder.Free();
       Modulator.Free();
       free(ModulatorOutput); ModulatorOutput=0;
       RateConverter.Free();
       free(ConverterOutput); ConverterOutput=0; }

   // preset internal arrays according to primary paramaters
   int Preset(MFSK_Parameters<Type> *NewParameters)
     {
       Parameters=NewParameters;

       BitsPerSymbol=Parameters->BitsPerSymbol;
       SymbolsPerBlock=Parameters->SymbolsPerBlock;

       // preset the input character buffer
       Input.Len=1024;
       if(Input.Preset()<0) goto Error;
       Monitor.Len=256;
       if(Monitor.Preset()<0) goto Error;

       // preset the encoder
       Encoder.BitsPerSymbol=BitsPerSymbol;
       if(Encoder.Preset()<0) goto Error;

       // preset the modulator
       if(Modulator.Preset(Parameters)<0) goto Error;

       if(ReallocArray(&ModulatorOutput,Modulator.OutputLen)<0) goto Error;

       // preset the rate converter
       RateConverter.OutputRate=Parameters->OutputSampleRate/Parameters->SampleRate;
       if(RateConverter.Preset()<0) goto Error;

       MaxOutputLen=(int)ceil(Parameters->SymbolSepar*Parameters->OutputSampleRate/Parameters->SampleRate+2);
       if(ReallocArray(&ConverterOutput,MaxOutputLen)<0) goto Error;

       Reset();

       return 0;

       Error: Free(); return -1; }

   void Reset(void)
     { Input.Reset();
       Monitor.Reset();
       SymbolPtr=0;
       State=0;
       RateConverter.Reset(); }

   // start the transmission
   void Start(void)
     { State|=State_Running; }

   // request to stop (and complete) the transmission
   // but the transmitter will only stop after transmitting all the data
   void Stop(void)
     { State|=State_StopReq; }

   // check if the transmission is still running (= not yet complete)
   int Running(void) const
     { return State&State_Running; }

   // put the character into the transmitter input queue
   int PutChar(uint8_t Char)
     { return Input.Write(Char); }

   // get one character from the monitor buffer
   int GetChar(uint8_t &Char)
     { return Monitor.Read(Char); }

   // get out the transmitter output (audio)
   int Output(Type *&OutputPtr)
     { if(SymbolPtr==0)                            // when at the block boundary
       { if((State&State_StopReq)&&Input.Empty()) // if the user requested to stop and no more characters
         { State=0; }                              // then simply stop
         else if(State&State_Running)             // otherwise when state is "running" then keep going
         { int Idx;                             // form a new block
           for(Idx=0; Idx<BitsPerSymbol; Idx++)   // loop over characters in a block
           { uint8_t Char;
             if(Input.Read(Char)<=0) break;       // get character from the input FIFO
             InputBlock[Idx]=Char;                 // put it into the block
             Monitor.Write(Char); }                // put this character into the monitor FIFO
           for(     ; Idx<BitsPerSymbol; Idx++)    // fill the unused part of the block
             InputBlock[Idx]=0;
           Encoder.EncodeBlock(InputBlock);         // encode the new block
         }
       }
       if(State&State_Running)                      // if state is "running" then
       { Modulator.Send(Encoder.OutputBlock[SymbolPtr]); // send out the next symbol of encoded block through the modulator
         SymbolPtr+=1; if(SymbolPtr>=SymbolsPerBlock) SymbolPtr=0; }
       int ModLen=Modulator.Output(ModulatorOutput);     // get the modulator output
       int ConvLen=RateConverter.Process(ModulatorOutput,ModLen,ConverterOutput); // correct the sampling rate
       if(ConvLen<0) return ConvLen;
       OutputPtr=ConverterOutput;
       return ConvLen; }

   int Output(int16_t *Buffer)
     { Type *OutputPtr;
       int OutputLen=Output(OutputPtr);
       ConvertToS16(ConverterOutput,Buffer,OutputLen);     // convert to 16-bit signed format (for the soundcard)
       return OutputLen; }

} ;

// =====================================================================

template <class Type=float>
 class MFSK_Synchronizer
{ public:

   MFSK_Parameters<Type> *Parameters;

  public:

  private:

   int FreqOffsets;                     // number of possible frequency offsets
   int BlockPhases;                     // number of possible time-phases within the FEC block
   MFSK_SoftDecoder<Type,Type> *Decoder;   // array of decoders
  public:
   int BlockPhase;                      // current running block time-phase
  private:
   CircularBuffer< LowPass3_Filter<Type>, int >  SyncNoiseEnergy; // FEC noise integrators
   CircularBuffer< LowPass3_Filter<Type>, int >  SyncSignal;      // FEC signal integrators
   Type SyncFilterWeight;                                    // weight for the integrators

  public:

   Type SyncBestSignal;                   // best signal
   int SyncBestBlockPhase;             // time-phase of the best signal        [FFT bins]
   int SyncBestFreqOffset;             // frequency offset of the best signal  [FFT spectral slices]
   Type SyncSNR;                          // S/N corresponding to the SyncBestSignal
   int DecodeReference;                   // when 0 then right in the middle of a FEC block

   Type PreciseFreqOffset;                // precise interpolated frequency offset [FFT bins]
   Type PreciseBlockPhase;                //                   and the FEC block phase [spectral slices]
   int StableLock;                     // is 1 when the sync. looks stable
   LowPass3_Filter<Type> FreqDrift;       // frequency drift rate [FFT bin / FEC block]
   LowPass3_Filter<Type> TimeDrift;       // block phase (time) drift rate [ppm]

  public:

   MFSK_Synchronizer()
     { Init(); }

   ~MFSK_Synchronizer()
     { Free(); }

   void Init(void)
     { Decoder=0; }

   void Free(void)
     { if(Decoder)
       { int Idx;
         for(Idx=0; Idx<FreqOffsets; Idx++)
           Decoder[Idx].Free();
         free(Decoder); Decoder=0;
       }
       SyncSignal.Free();
       SyncNoiseEnergy.Free();
     }

   // resize internal arrays according the parameters
   int Preset(MFSK_Parameters<Type> *NewParameters)
     { Parameters=NewParameters;
	 
	   int Idx;

       if(Decoder)
       { for(Idx=0; Idx<FreqOffsets; Idx++)
           Decoder[Idx].Free();
       }

       FreqOffsets=2*Parameters->RxSyncMargin*Parameters->CarrierSepar+1;
       BlockPhases=Parameters->SpectraPerSymbol*Parameters->SymbolsPerBlock;

       if(ReallocArray(&Decoder,FreqOffsets)<0) goto Error;
       for(Idx=0; Idx<FreqOffsets; Idx++)
         Decoder[Idx].Init();
       for(Idx=0; Idx<FreqOffsets; Idx++)
         if(Decoder[Idx].Preset(Parameters)<0) goto Error;

       SyncSignal.Width=FreqOffsets;
       SyncSignal.Len=BlockPhases;
       if(SyncSignal.Preset()<0) goto Error;

       SyncNoiseEnergy.Width=FreqOffsets;
       SyncNoiseEnergy.Len=BlockPhases;
       if(SyncNoiseEnergy.Preset()<0) goto Error;

       SyncFilterWeight = 1.0/Parameters->RxSyncIntegLen;

       Reset();

       return 0;

       Error: Free(); return -1; }

   void Reset(void)
     { int Idx;

       for(Idx=0; Idx<FreqOffsets; Idx++)
         Decoder[Idx].Reset();

       SyncSignal.Clear();
       SyncNoiseEnergy.Clear();

       BlockPhase=0;

       SyncBestSignal=0;
       SyncBestBlockPhase=0;
       SyncBestFreqOffset=0;
       SyncSNR=0;
       DecodeReference=(-BlockPhases/2);

       PreciseFreqOffset=0;
	   PreciseBlockPhase=0;
       StableLock=0;
       FreqDrift=0;
       TimeDrift=0;

	 }

   void Process(Type *Spectra)
     { int Offset;

	   MFSK_SoftDecoder<Type,Type> *DecoderPtr = Decoder;
       LowPass3_Filter<Type> *SignalPtr        = SyncSignal[BlockPhase];
       LowPass3_Filter<Type> *NoiseEnergyPtr   = SyncNoiseEnergy[BlockPhase];

       // printf("%3d:",BlockPhase);
       Type BestSliceSignal=0;
       int BestSliceOffset=0;
	   for(Offset=0; Offset<FreqOffsets; Offset++)
	   { DecoderPtr->SpectralInput(Spectra+Offset);
	     DecoderPtr->Process();
         Type NoiseEnergy = DecoderPtr->NoiseEnergy;
         Type Signal = DecoderPtr->Signal;

         // printf(" %4.1f",Signal);

         NoiseEnergyPtr->Process(NoiseEnergy, SyncFilterWeight);
         SignalPtr->Process(Signal, SyncFilterWeight);
         Signal=SignalPtr->Output;

         // printf("/%4.1f",Signal);

         if(Signal>BestSliceSignal)
         { BestSliceSignal=Signal;
           BestSliceOffset=Offset; }

         DecoderPtr++;
         NoiseEnergyPtr++;
         SignalPtr++;
	   } // printf("\n");

       if(BlockPhase==SyncBestBlockPhase)
       { SyncBestSignal=BestSliceSignal;
         SyncBestFreqOffset=BestSliceOffset;
       }
       else
       { if(BestSliceSignal>SyncBestSignal)
         { SyncBestSignal=BestSliceSignal;
           SyncBestBlockPhase=BlockPhase;
           SyncBestFreqOffset=BestSliceOffset; }
       }
/*
       printf("MFSK_Synchronizer: %4.1f @ %3d:%3d\n",
	                             SyncBestSignal, SyncBestBlockPhase, SyncBestFreqOffset);
*/
       DecodeReference=(int)BlockPhase-(int)SyncBestBlockPhase;
       if(DecodeReference<0) DecodeReference+=BlockPhases;
	   DecodeReference-=(int)(BlockPhases/2);
       if(DecodeReference==0)
       { Type BestNoise=SyncNoiseEnergy[SyncBestBlockPhase][SyncBestFreqOffset].Output;
         if(BestNoise>0) BestNoise=sqrt(BestNoise);
		            else BestNoise=0;
         const Type MinNoise=(Type)Parameters->SymbolsPerBlock/10000;
         if(BestNoise<MinNoise) BestNoise=MinNoise;

         SyncSNR = SyncBestSignal/BestNoise;

         Type NewPreciseFreqOffset;
         Type SignalPeak;
         LowPass3_Filter<Type> *Signal = SyncSignal[SyncBestBlockPhase];
         int FitIdx=Limit(SyncBestFreqOffset,(int)1,(int)(FreqOffsets-2));
		 int FitOK=FitPeak(NewPreciseFreqOffset, SignalPeak,
		                   Signal[FitIdx-1].Output, Signal[FitIdx].Output, Signal[FitIdx+1].Output);
         if(FitOK<0) NewPreciseFreqOffset=SyncBestFreqOffset;
                else NewPreciseFreqOffset=FitIdx+Limit(NewPreciseFreqOffset,(Type)-1.0,(Type)1.0);

         Type NewPreciseBlockPhase;
         int FitIdxL=SyncBestBlockPhase;
	     SyncSignal.DecrPtr(FitIdxL);
         int FitIdxC=SyncBestBlockPhase;
         int FitIdxR=SyncBestBlockPhase;
	     SyncSignal.IncrPtr(FitIdxR);
		 FitOK=FitPeak(NewPreciseBlockPhase, SignalPeak,
		                 SyncSignal[FitIdxL][SyncBestFreqOffset].Output,
		                 SyncSignal[FitIdxC][SyncBestFreqOffset].Output,
		                 SyncSignal[FitIdxR][SyncBestFreqOffset].Output);
         if(FitOK<0) { NewPreciseBlockPhase=SyncBestBlockPhase; }
                else { NewPreciseBlockPhase+=FitIdxC; SyncSignal.WrapPhase(NewPreciseBlockPhase); }

         Type FreqDelta=NewPreciseFreqOffset-PreciseFreqOffset;
         Type PhaseDelta=NewPreciseBlockPhase-PreciseBlockPhase;
         SyncSignal.WrapDiffPhase(PhaseDelta);

         Type DeltaDist2=FreqDelta*FreqDelta+PhaseDelta*PhaseDelta;
         if((DeltaDist2<=1.0)&&(SyncSNR>=Parameters->RxSyncThreshold))
		 { StableLock=1;
           FreqDrift.Process(FreqDelta, SyncFilterWeight);
           TimeDrift.Process(PhaseDelta/BlockPhases, SyncFilterWeight); }
		 else
		 { StableLock=0; FreqDrift=0; TimeDrift=0; }
/*
         printf("%d: %4d (%6.2f %+5.2f %+5.0f ppm) / %+2d (%+5.3f %+5.2f %+7.4f) => %4.1f/%3.1f = %4.1f\n",
             StableLock,
             SyncBestBlockPhase, NewPreciseBlockPhase, PhaseDelta, 1E6*TimeDrift.Output,
             (int)SyncBestFreqOffset-(FreqOffsets/2),  NewPreciseFreqOffset-(FreqOffsets/2), FreqDelta, FreqDrift.Output,
             SyncBestSignal, BestNoise,
             SyncSNR);
*/
         PreciseFreqOffset=NewPreciseFreqOffset;
         PreciseBlockPhase=NewPreciseBlockPhase;

       }
       SyncSignal.IncrPtr(BlockPhase);
    }

   Type FEC_SNR(void)
     { return SyncSNR; }

   Type FrequencyOffset(void)
     { return (PreciseFreqOffset-(FreqOffsets/2))*Parameters->FFTbinBandwidth(); }

   Type FrequencyDriftRate(void)
     { return FreqDrift.Output*Parameters->FFTbinBandwidth()/Parameters->BlockPeriod(); }

   Type TimeDriftRate(void)
     { return TimeDrift.Output; }

  } ;

// =====================================================================

/*

How to use the MFSK_Receiver class:

1. create an object like:

   #include "mfsk.h"

   MFSK_Receiver<float> Receiver;

2. Set the parameters, for example:

   Receiver.Tones           = 32;     // number of tones (symbols)
   Receiver.Bandwidth       = 1000;   // bandwidth [Hz]
   Receiver.SyncMargin      = 8;      // synchronizer tune margin [tone freq. spacing]
   Receiver.SyncIntegLen    = 4;      // synchronizer integration period [FEC blocks]
   Receiver.SyncThreshold   = 3.2;    // S/N threshold for printing
   Receiver.SampleRate      = 8000.0; // internal processor sampling rate [Hz]
   Receiver.InputSampleRate = 8000.0; // soundcard sampling rate [Hz]

   You don't need to set all the parameters, as upon creation
   of the Receiver object they are already given certain default
   values.

   If you changed parameters at one time and want later to go back
   to the default values you can call: Receiver.Default();

3. Preset the Receiver internal arrays for the parameters you just set:

   if(Receiver.Preset()<0)
     printf("Not enough RAM or another problem\n");

   Each time you change the parameters you need to call Preset()
   in order to resize the internal arrays. Preset() will as well
   destroy all data being in the process of decoding, if you need
   this data then call first Receiver.Flush()

4. Read back the parameters you set in point 1., they could have been adjusted
   by Preset() to their closest allowed values.

5. Feed the audio into the Receiver:

   Receiver.Process(AudioBuffer, BufferLength);

   AudioBuffer can be an array of int16_t (16-bit signed integers)
   that you fill with the data from the soundcard. I suggest you feed
   the receiver with batches of 512 or 1024 samples, but in can be any number
   of samples at a time.

6. Call GetChar(Char) to get decoded characters. Note, that
   characters come in batches, and so, you need to call GetChar()
   possibly several times. GetChar() returns 0 when the decoder FIFO
   is empty or 1 when not empty. In the latter case the argument
   contains the character read form the FIFO. The loop could be like:
   
   for( ; ; )
   { uint8_t Char;
     if(Receiver.GetChar(Char)==0) break;
     printf("%c",Char);
   }

   Keep in mind that you may see (random) control code characters here,
   and thus you must be able to deal with them. I suggest to process
   only carriage return (code=13) and Backspace (code=8). NUL (code=0)
   is the idle character: it is being sent when there is no text to be sent.

7. At any time you can read the signal-to-noise ratio of the
   incoming signal by calling Receiver.SignalToNoiseRatio() or frequency offset
   by calling Receiver.FrequencyOffset()

8. When the user decides to turn off the receiver and switch over
   to transmitt you may still call Receiver.Flush()in order to flush
   the data still being buffered in the decoder pipeplines.


*/

template <class Type=float>
 class MFSK_Receiver
{ public:

   MFSK_Parameters<Type> *Parameters;

  private:

   RateConverter<Type> RateConverter;
   Seq<Type> InputBuffer;
   MFSK_InputProcessor<Type> InputProcessor; // equalizes the input spectrum
                                             // and removes coherent interferences
   MFSK_Demodulator<Type> Demodulator;       // spectral (FFT) demodulator
   MFSK_Synchronizer<Type> Synchronizer;     // synchronizer
   MFSK_SoftIterDecoder<Type> Decoder;       // iterative decoder
   FIFO<uint8_t, int> Output;                     // buffer for decoded characters

  public:

   MFSK_Receiver()
     { Init(); }

   ~MFSK_Receiver()
     { Free(); }

   void Init(void)
     { }

   void Free(void)
     { RateConverter.Free();
       InputBuffer.Free();
       InputProcessor.Free();
       Demodulator.Free();
       Synchronizer.Free();
       Decoder.Free();
       Output.Free(); }

   // resize internal arrays according the parameters
   int Preset(MFSK_Parameters<Type> *NewParameters)
     { Parameters=NewParameters;

       RateConverter.OutputRate=Parameters->SampleRate/Parameters->InputSampleRate;
       if(RateConverter.Preset()<0) goto Error;

       InputProcessor.WindowLen=32*Parameters->SymbolSepar;
       if(InputProcessor.Preset()<0) goto Error;

       if(InputBuffer.EnsureSpace(InputProcessor.WindowLen+2048)<0) goto Error;

       if(Demodulator.Preset(Parameters)<0) goto Error;
       if(Synchronizer.Preset(Parameters)<0) goto Error;
       if(Decoder.Preset(Parameters)<0) goto Error;

       Output.Len=1024;
       if(Output.Preset()<0) goto Error;

       return 0;

       Error: Free(); return -1; }

   void Reset(void)
     { RateConverter.Reset();
       InputBuffer.Clear();
       InputProcessor.Reset();
       Demodulator.Reset();
       Synchronizer.Reset();
       Output.Reset(); }

   Type SyncSNR(void)
   { return Synchronizer.FEC_SNR(); }

   Type FrequencyOffset(void)
   { return Synchronizer.FrequencyOffset(); }

   Type FrequencyDrift(void)
   { return Synchronizer.FrequencyDriftRate(); }

   Type TimeDrift(void)
   { return Synchronizer.TimeDriftRate(); }

   Type InputSNRdB(void)
   { return Decoder.InputSNRdB(); }

   // process an audio batch: first the input processor, then the demodulator
   template <class InpType>
    int Process(InpType *Input, int InputLen)
     { if(RateConverter.Process(Input, InputLen, InputBuffer)<0) return -1;
       ProcessInputBuffer();
	   return 0; }

   void Flush(void)
     { ProcessInputBuffer();
       int Idx;

	   for(Idx=InputBuffer.Len; Idx<InputProcessor.WindowLen; Idx++)
	     InputBuffer[Idx]=0;
	   InputBuffer.Len=InputProcessor.WindowLen;
       ProcessInputBuffer();

	   for(Idx=0; Idx<InputProcessor.WindowLen; Idx++)
	     InputBuffer[Idx]=0;
       int FlushLen=Parameters->SymbolSepar*Parameters->SymbolsPerBlock*Parameters->RxSyncIntegLen*2;
       for(Idx=0; Idx<FlushLen; Idx+=InputProcessor.WindowLen)
       { InputBuffer.Len=InputProcessor.WindowLen;
         ProcessInputBuffer(); }
     }

   // get one character from the output buffer
   int GetChar(uint8_t &Char)
     { return Output.Read(Char); }

  private:

   // process the input buffer: first the input processor, then the demodulator
    void ProcessInputBuffer(void)
     { while(InputBuffer.Len>=InputProcessor.WindowLen)
       { InputProcessor.Process(InputBuffer.Elem);
         InputBuffer.Delete(0,InputProcessor.WindowLen);
         int Idx;
         for(Idx=0; Idx<InputProcessor.WindowLen; Idx+=Parameters->SymbolSepar )
           ProcessSymbol(InputProcessor.Output+Idx);
       }
     }

   // process (through the demodulator) an audio batch corresponding to one symbol
   // (demodulator always works with audio batches corresponding to one symbol period)
   template <class InpType>
    void ProcessSymbol(InpType *Input)
   { int SpectraPerSymbol=Parameters->SpectraPerSymbol;
     Demodulator.Process(Input);
    int HistOfs;
    for(HistOfs=(-SpectraPerSymbol); HistOfs<0; HistOfs++)
    { Type *Spectra = Demodulator.HistoryPtr(HistOfs);
      Synchronizer.Process(Spectra);
      int SpectraPerBlock=Parameters->SpectraPerBlock;
      if(Synchronizer.DecodeReference==0)
      { 
/*
	    printf("%d: SNR=%4.1f, %+5.1f Hz, %d-%d\n",
	           Synchronizer.StableLock, SyncSNR(), FrequencyOffset(),
			   Synchronizer.BlockPhase, Synchronizer.SyncBestBlockPhase);
*/
        if(Synchronizer.StableLock)
        { int TimeOffset = (HistOfs-((Parameters->RxSyncIntegLen+1)*SpectraPerBlock+SpectraPerBlock/2-1));
          int FreqOffset = Synchronizer.SyncBestFreqOffset;

          Type BestSignal=0;
          int BestTime=0;
          int BestFreq=0;
          int FreqSearch;
          for(FreqSearch=(-1); FreqSearch<=1; FreqSearch++)
          { int TimeSearch;
            for(TimeSearch=(-2); TimeSearch<=2; TimeSearch++)
            { int Error=Demodulator.PickBlock(Decoder.Input, TimeOffset+TimeSearch,FreqOffset+FreqSearch);
              if(Error<0) continue;
	      Decoder.Process(8);
              // printf("%+2d/%+2d: ", FreqSearch, TimeSearch);
              // Decoder.PrintSNR();
              Type Signal=Decoder.Input_SignalEnergy;
              if(Signal>BestSignal)
              { BestSignal=Signal; BestFreq=FreqSearch; BestTime=TimeSearch; }
            }
          }

          Demodulator.PickBlock(Decoder.Input, TimeOffset+BestTime,FreqOffset+BestFreq);
          Decoder.Process(32);
          // printf("Best: %+2d/%+2d: ", BestFreq, BestTime);
          // Decoder.PrintSNR();
          Decoder.WriteOutputBlock(Output);

        }

      }
    }

   }

} ;

// =====================================================================

#endif // of __MFSK_H__
