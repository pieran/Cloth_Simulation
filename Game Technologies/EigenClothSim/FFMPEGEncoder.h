//==============================================================================//
// INDUSTRIAL MONITORING AND CONTROL 2014                                       //
//==============================================================================//
//      ___           ___           ___                                         // 
//     /\__\         /\  \         /\__\                                        //
//    /:/  /        |::\  \       /:/  /                                        //
//   /:/__/         |:|:\  \     /:/  /                                         //
//  /::\  \ ___   __|:|\:\  \   /:/  /  ___                                     //
// /:/\:\  /\__\ /::::|_\:\__\ /:/__/  /\__\                                    //
// \/__\:\/:/  / \:\~~\  \/__/ \:\  \ /:/  /                                    //
//      \::/  /   \:\  \        \:\  /:/  /                                     //
//      /:/  /     \:\  \        \:\/:/  /                                      //
//     /:/  /       \:\__\        \::/  /                                       //
//     \/__/         \/__/         \/__/                                        //
//                                                                              //
// Webstore and Free Code:                                                      //
// http://www.imc-store.com.au/                                                 //
//                                                                              //
// Industrial/Commercial Information                                            //
// http://imcontrol.com.au/                                                     //
//                                                                              //
//------------------------------------------------------------------------------//
// FFMPEG ENCODER CLASS EXAMPLE                                                 //
//------------------------------------------------------------------------------//
// Notes: This C++ program is designed to simplify talking with the FFMPEG      //
// library. It provides some basic functions to encode a video from your code.  //
//                                                                              //
// Usage: Use this project as a template for your own program, or add the       //
// FFMPEG class to your project                                                 //
//==============================================================================//


#pragma once
#define _CRT_SECURE_NO_WARNINGS

//=============================
// Includes
//-----------------------------
// FFMPEG is writen in C so we need to use extern "C"
//-----------------------------
extern "C" {
	#define INT64_C(x) (x ## LL)
	#define UINT64_C(x) (x ## ULL)

	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#include <math.h>
	#include <libavutil/opt.h>
	#include <libavutil/mathematics.h>
	#include <libavformat/avformat.h>
	#include <libswscale/swscale.h>
	#include <libswresample/swresample.h>
	#include <libavutil/imgutils.h>
	#include <libavcodec/avcodec.h>

	static int sws_flags = SWS_BICUBIC;
}

//=============================
// FFMPEG Class
//-----------------------------
class FFMPEGEncoder
{
public:
	//-----------------------------
	//Define our subroutines that we want the outside world to access
	//-----------------------------

	FFMPEGEncoder();
	~FFMPEGEncoder();

	void SetupVideo(const char * filename, int Width, int Height, int FPS, int GOB, int BitPerSecond);
	void WriteDummyFrame();
	void WriteFrame(const char * RGBFrame);
	void CloseVideo(void);

	int  GetVideoWidth(void)		{ return m_AVIMOV_WIDTH; }
	int  GetVideoHeight(void)	{ return m_AVIMOV_HEIGHT; }
	bool GetIsRecording()		{ return m_AVIMutex; }
private:
	//-----------------------------
	//Define our subroutines that only our class can access
	//-----------------------------

	int		m_sws_flags;
	int    m_AVIMOV_FPS;
	int    m_AVIMOV_GOB;
	int    m_AVIMOV_BPS;
	int		m_frame_count;
	int    m_AVIMOV_WIDTH;
	int    m_AVIMOV_HEIGHT;

	bool m_AVIMutex;

	double m_video_time;

	char *m_filename;

	void CloseCodec(void);
	void WriteFrame(void);
	void SetupCodec(const char *filename, int codec_id);

	AVFrame *m_frame;
	AVCodecContext *m_c;
	AVStream *m_video_st;
	AVOutputFormat *m_fmt;
	AVFormatContext *m_oc;
	AVCodec *m_video_codec;
	AVPicture m_src_picture, m_dst_picture;
};