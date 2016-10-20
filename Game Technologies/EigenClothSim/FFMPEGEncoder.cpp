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

#pragma warning(disable:4996)

//=============================
// Includes
//-----------------------------
// FFMPEG is writen in C so we need to use extern "C"
//-----------------------------
#include "FFMPEGEncoder.h"



//=============================
// Our Creator
//-----------------------------
// Create any necessary blocks of memory and init variables
//-----------------------------
FFMPEGEncoder::FFMPEGEncoder()
{
	m_filename = new char[1024];
	m_AVIMutex = false;
}



//=============================
// Our Destructor
//-----------------------------
// Delete any allocated memory
//-----------------------------
FFMPEGEncoder::~FFMPEGEncoder()
{
	delete[] m_filename;
}




//=============================
// Setup Video Codec
//-----------------------------
// Sets up the video codec
//-----------------------------
void FFMPEGEncoder::SetupCodec(const char *filename, int codec_id)
{
	//---------------------
	//The basic flow for opening a video is:
	//---------------------
	//Register FFMPEG -> Allocate Media Context -> Allocate Stream Format Context ->
	//Allocate Codec Context -> Allocate Frames -> Open Media File -> Start Writing Frames!!

	//If file is already open then don't allow user to open it again
	if (m_AVIMutex) { return; }

	//Some inits
	int ret;
	m_sws_flags = SWS_BICUBIC;
	m_frame_count = 0;

	//You must call these subroutines otherwise you get errors!!
	avcodec_register_all();
	av_register_all();

	//allocate the output media context
	avformat_alloc_output_context2(&m_oc, NULL, NULL, filename);
	if (!m_oc) {
		//if context failed check by specifying container
		avformat_alloc_output_context2(&m_oc, NULL, "avi", filename);
	}
	if (!m_oc) {
		//if context failed check by specifying container a bit further
		avformat_alloc_output_context2(&m_oc, NULL, "avi", "c:\\temp.avi");
	}
	if (!m_oc) {
		//If you reach this point something's wrong
		printf("ERROR\n");
		return;
	}

	//Get the format determined by the conatiner
	m_fmt = m_oc->oformat;

	// Add the audio and video streams using the default format codecs
	// and initialize the codecs.
	m_video_st = NULL;

	//Setup the codecs
	m_fmt->video_codec = AV_CODEC_ID_H264;
	m_fmt->audio_codec = AV_CODEC_ID_NONE;

	// Add an output stream.
	{
		//Init AV Stream
		AVStream *st;

		//find the encoder
		m_video_codec = avcodec_find_encoder(m_fmt->video_codec);
		if (!(m_video_codec)) {
			printf("ERROR\n");
			return;
		}

		//Create new video stream
		st = avformat_new_stream(m_oc, m_video_codec);
		if (!st) {
			printf("ERROR\n");
			return;
		}
		st->id = m_oc->nb_streams - 1;

		st->time_base.den = m_AVIMOV_FPS;
		st->time_base.num = 1;

		//Set codec context
		m_c = st->codec;

		//Setup fundumental video stream parameters
		m_c->codec_id = m_fmt->video_codec;
		m_c->bit_rate = m_AVIMOV_BPS;                   //Bits Per Second
		m_c->width = m_AVIMOV_WIDTH;                 //Note Resolution must be a multiple of 2!!
		m_c->height = m_AVIMOV_HEIGHT;         //Note Resolution must be a multiple of 2!!
		//m_c->time_base.den = m_AVIMOV_FPS;       //Frames per second
		//m_c->time_base.num = 1;
		m_c->gop_size = m_AVIMOV_GOB;       // Intra frames per x P frames
		m_c->pix_fmt = AV_PIX_FMT_YUV420P;//Do not change this, H264 needs YUV format not RGB


										  //Some formats want stream headers to be separate.
		if (m_oc->oformat->flags & AVFMT_GLOBALHEADER)
			m_c->flags |= CODEC_FLAG_GLOBAL_HEADER;

		//Set our video stream pointer
		m_video_st = st;
	}


	// Now that all the parameters are set, we can open the audio and
	// video codecs and allocate the necessary encode buffers.
	{
		//Allocated Codec Context
		AVCodecContext *c = m_video_st->codec;

		//Open the codec
		ret = avcodec_open2(c, m_video_codec, NULL);
		if (ret < 0) {
			printf("ERROR\n");
			return;
		}

		//allocate and init a re-usable frame
		m_frame = av_frame_alloc();
		if (!m_frame) {
			printf("ERROR\n");
			return;
		}

		//Allocate the encoded raw picture.
		ret = avpicture_alloc(&m_dst_picture, c->pix_fmt, c->width, c->height);
		if (ret < 0) {
			printf("ERROR\n");
			return;
		}

		//Allocate RGB frame that we can pass to the YUV frame
		ret = avpicture_alloc(&m_src_picture, AV_PIX_FMT_RGB24, c->width, c->height);
		if (ret < 0) {
			printf("ERROR\n");
			return;
		}

		//Copy data and linesize picture pointers to frame
		*((AVPicture *)m_frame) = m_dst_picture;
		m_frame->format = c->pix_fmt;
		m_frame->width = c->width;
		m_frame->height = c->height;
	}

	//Tell FFMPEG that we are going to write encoded frames to a file
	av_dump_format(m_oc, 0, filename, 1);

	//open the output file, if needed
	if (!(m_fmt->flags & AVFMT_NOFILE)) {
		ret = avio_open(&m_oc->pb, filename, AVIO_FLAG_WRITE);
		if (ret < 0) {
			printf("ERROR\n");
			return;
		}
	}


	//Write the stream header, if any.
	ret = avformat_write_header(m_oc, NULL);


	if (ret < 0) {
		printf("ERROR\n");
		return;
	}

	//Set frame count to zero
	if (m_frame)
		m_frame->pts = 0;

	//Frame is initalised!
	m_AVIMutex = true;

}

//=============================
// Write Video Frame
//-----------------------------
// Writes a Video Frame that is stored in m_src_picture
//-----------------------------
void FFMPEGEncoder::WriteFrame(void) {

	//If video is not initalised then don't write frame
	if (!m_AVIMutex) { return; }

	//Calculate video time    
	m_video_time = m_video_st ? m_video_st->pts.val * av_q2d(m_video_st->time_base) : 0.0;

	//Inits
	int ret;
	static struct SwsContext *sws_ctx;
	AVCodecContext *c = m_video_st->codec;

	//If we haven't already done so init the context of the frame conversion from RGB to YUV
	if (!sws_ctx) {
		sws_ctx = sws_getContext(c->width, c->height, AV_PIX_FMT_RGB24,
			c->width, c->height, AV_PIX_FMT_YUV420P,
			sws_flags, NULL, NULL, NULL);
		if (!sws_ctx) {
			printf("ERROR\n");
			return;
		}
	}

	//Convert RGB frame (m_src_picture) to and YUV frame (m_dst_picture)
	sws_scale(sws_ctx,
		m_src_picture.data, m_src_picture.linesize,
		0, c->height, m_dst_picture.data, m_dst_picture.linesize);


	//Some inits for encoding the frame
	AVPacket pkt = { 0 };
	int got_packet;
	av_init_packet(&pkt);

	//Encode the frame
	ret = avcodec_encode_video2(c, &pkt, m_frame, &got_packet);
	if (ret < 0) {
		printf("ERROR\n");
		return;
	}

	//If size of encoded frame is zero, it means the image was buffered.
	if (!ret && got_packet && pkt.size) {
		pkt.stream_index = m_video_st->index;
		
		//Write the compressed frame to the media file.
		ret = av_interleaved_write_frame(m_oc, &pkt);

		//if non-zero then it means that there was something wrong writing the frame to
		//the file
		if (ret != 0) {
			printf("ERROR\n");
			return;
		}
	}
	else {
		ret = 0;
	}

	//Increment Frame counter
	m_frame_count++;
	m_frame->pts += 1;
}



//=============================
// Setup Video
//-----------------------------
// Sets up the Video and Opens the Video File
//-----------------------------
void FFMPEGEncoder::SetupVideo(const char * filename, int Width, int Height, int FPS, int GOB, int BitPerSecond) {
	//Copy filename to local string
	sprintf(m_filename, filename);

	//Set movie parameters
	m_AVIMOV_WIDTH = Width;      //Movie width
	m_AVIMOV_HEIGHT = Height;    //Movie height
	m_AVIMOV_FPS = FPS;          //Movie frames per second
	m_AVIMOV_GOB = GOB;          //I frames per no of P frames, see note below!
	m_AVIMOV_BPS = BitPerSecond; //Bits per second, if this is too low then movie will become garbled

								 //Pass Parameters to Setup Codec as H264 file
	SetupCodec(m_filename, AV_CODEC_ID_H264);

	//------------
	//NOTE: on GOB
	//------------
	//An I frame is an entire frame stored. A P frame is only a partial frame that is
	//encoded based on the previous frame.
	//
	//Think of I frames as refreshing the entire video, and P frames as just encoding the
	//moving bits of the frame.
	//
	//A GOB of 10 means that for every 10 P frames there is only 1 I frame, i.e. the frame is
	//only 'refreshed' every 10 frames.
	//------------
}

//=============================
// Close Codec
//-----------------------------
// Close Video Codec
//-----------------------------
void FFMPEGEncoder::CloseCodec(void) {

	//If video is not initalised then don't close frame
	if (!m_AVIMutex) { return; }

	//Write trailing bits
	av_write_trailer(m_oc);

	//Close Video codec
	avcodec_close(m_video_st->codec);

	//Free our frames
	av_free(m_src_picture.data[0]);
	av_free(m_dst_picture.data[0]);
	av_free(m_frame);

	//A bit more cleaning
	if (!(m_fmt->flags & AVFMT_NOFILE))
		avio_close(m_oc->pb);

	avformat_free_context(m_oc);

	//Set open flag clear
	m_AVIMutex = false;
}


//=============================
// Close Codec
//-----------------------------
// Close Video Codec
//-----------------------------
void FFMPEGEncoder::CloseVideo(void) {

	//Close Video Codec and write end of file
	CloseCodec();
}



//=============================
// Write Dummy Frame
//-----------------------------
// Processes a Randomly Generated RGB frame
//-----------------------------
void FFMPEGEncoder::WriteDummyFrame() {

	//Step through height of frame
	for (int y = 0; y<m_c->height; y++) {  //Height Loop

										   //Step through width of frame
		for (int x = 0; x<m_c->width; x++) { //Width Loop

											 //Save RGB frame to FFMPEG's source frame
			m_src_picture.data[0][y * m_src_picture.linesize[0] + x * 3 + 0] = char(rand() * 255);  //Red Channel
			m_src_picture.data[0][y * m_src_picture.linesize[0] + x * 3 + 1] = char(rand() * 255);  //Green Channel
			m_src_picture.data[0][y * m_src_picture.linesize[0] + x * 3 + 2] = char(rand() * 255);  //Blue Channel
		}
	}

	//Send frame off to FFMPEG for encoding
	WriteFrame();
}



//=============================
// Write Frame
//-----------------------------
// Processes an RGB frame supplied by the user
//-----------------------------
void FFMPEGEncoder::WriteFrame(const char * RGBFrame) {

	//Data should be in the format RGBRGBRGB...etc and should be Width*Height*3 long

	//Step through height of frame
	for (int y = 0; y<m_c->height; y++) {  //Height Loop

										   //Step through width of frame
		for (int x = 0; x<m_c->width; x++) { //Width Loop

											 //Save RGB frame to FFMPEG's source frame
			m_src_picture.data[0][y * m_src_picture.linesize[0] + x * 3 + 0] = RGBFrame[(y*m_c->width + x) * 3 + 0];  //Red Channel
			m_src_picture.data[0][y * m_src_picture.linesize[0] + x * 3 + 1] = RGBFrame[(y*m_c->width + x) * 3 + 1];  //Green Channel
			m_src_picture.data[0][y * m_src_picture.linesize[0] + x * 3 + 2] = RGBFrame[(y*m_c->width + x) * 3 + 2];  //Blue Channel
		}
	}

	//Send frame off to FFMPEG for encoding
	WriteFrame();
}