#include "VideoEncoder.h"

VideoEncoder::VideoEncoder()
{
	glGenBuffers(1, &m_CopyPixelBuffer);


	/* register all the codecs */
	av_register_all();
	avcodec_register_all();

	//av_dict_set(&opt, "qp", "0", 0);

	fmt = NULL;
}

VideoEncoder::~VideoEncoder()
{
	if (fmt) EndEncoding();

	glDeleteBuffers(1, &m_CopyPixelBuffer);
}

void VideoEncoder::BeginEncoding(int width, int height, const std::string& filename)
{
	if (fmt) EndEncoding();

	crop_x = width & 0x1;
	crop_y = height & 0x1;
	width -= crop_x;
	height -= crop_y;

	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, m_CopyPixelBuffer);
	glBufferDataARB(GL_PIXEL_PACK_BUFFER_ARB, width * height * 3, NULL, GL_STATIC_READ_ARB);
	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, 0);

	printf("Encode video file %s\n", filename.c_str());

	/* allocate the output media context */
	avformat_alloc_output_context2(&oc, NULL, NULL, filename.c_str());
	if (!oc) {
		printf("Could not deduce output format from file extension: using MPEG.\n");
		avformat_alloc_output_context2(&oc, NULL, "mpeg", filename.c_str());
	}
	if (!oc)
	{
		printf("Error: unable to create video stream!\n");
		system("pause");
		exit(1);
	}

	fmt = oc->oformat;

	/* Add the audio and video streams using the default format codecs
	* and initialize the codecs. */
	if (fmt->video_codec == AV_CODEC_ID_NONE)
	{
		printf("Error: unable to create video stream!\n");
		system("pause");
		exit(1);
	}

	add_stream(&video_st, oc, &video_codec, fmt->video_codec, width, height);

	/* Now that all the parameters are set, we can open the audio and
	* video codecs and allocate the necessary encode buffers. */
	open_video(oc, video_codec, &video_st, opt);

	av_dump_format(oc, 0, filename.c_str(), 1);

	/* open the output file, if needed */
	if (!(fmt->flags & AVFMT_NOFILE)) {
		int ret = avio_open(&oc->pb, filename.c_str(), AVIO_FLAG_WRITE);
		if (ret < 0) {
			fprintf(stderr, "Could not open '%s'\n", filename.c_str());
			system("pause");
			exit(1);
		}
	}

	/* Write the stream header, if any. */
	int ret = avformat_write_header(oc, &opt);
	if (ret < 0) {
		fprintf(stderr, "Error occurred when opening output file\n");
		system("pause");
		exit(1);
	}
}

void VideoEncoder::EndEncoding()
{
	if (fmt)
	{
		/* Write the trailer, if any. The trailer must be written before you
		* close the CodecContexts open when you wrote the header; otherwise
		* av_write_trailer() may try to use memory that was freed on
		* av_codec_close(). */
		av_write_trailer(oc);

		/* Close each codec. */
		close_stream(oc, &video_st);

		if (!(fmt->flags & AVFMT_NOFILE))
			/* Close the output file. */
			avio_closep(&oc->pb);

		/* free the stream */
		avformat_free_context(oc);
	}
	printf("Video Finished Encoding\n--------------------------\n\n");
	fmt = NULL;
}

void VideoEncoder::EncodeFrame()
{
	if (!fmt)
		return;


	//fflush(stdout);



	/*glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, m_CopyPixelBuffer);

	// set the target framebuffer to read
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, c->width, c->height, GL_RGB, GL_UNSIGNED_BYTE, 0);

	// map the PBO to process its data by CPU
	//glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, m_ReadPixelBuffer);
	GLubyte* ptr = (GLubyte*)glMapBufferARB(GL_PIXEL_PACK_BUFFER_ARB, GL_READ_ONLY_ARB);
	if (ptr)
	{
	RGB2Yuv420p(frame, ptr, c->width, c->height);
	glUnmapBufferARB(GL_PIXEL_PACK_BUFFER_ARB);
	}
	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, 0);*/




	write_video_frame(oc, &video_st);

}

void VideoEncoder::RGB2Yuv420p(AVFrame *frame, uint8_t *rgb, const int &width, const int &height)
{
	const size_t image_size = width * height;
	uint8_t *dst_y = frame->data[0];
	uint8_t *dst_u = frame->data[1];
	uint8_t *dst_v = frame->data[2];

	// Y plane
	for (size_t y = 0; y < height; y++) {
		for (size_t x = 0; x < width; x++) {
			const size_t i = y*frame->linesize[0] + x;
			const size_t j = (height - y - 1)*width + x;
			dst_y[i] = ((66 * rgb[3 * j] + 129 * rgb[3 * j + 1] + 25 * rgb[3 * j + 2]) >> 8) + 16;
		}
	}

	// U plane
	for (size_t y = 0; y < height; y += 2) {
		for (size_t x = 0; x < width; x += 2) {
			const size_t i = y / 2 * frame->linesize[1] + x / 2;
			const size_t j = (height - y - 1)*width + x;
			dst_u[i] = ((-38 * rgb[3 * j] + -74 * rgb[3 * j + 1] + 112 * rgb[3 * j + 2]) >> 8) + 128;
		}
	}
	// V plane
	for (size_t y = 0; y < height; y += 2) {
		for (size_t x = 0; x < width; x += 2) {
			const size_t i = y / 2 * frame->linesize[2] + x / 2;
			const size_t j = (height - y - 1)*width + x;
			dst_v[i] = ((112 * rgb[3 * j] + -94 * rgb[3 * j + 1] + -18 * rgb[3 * j + 2]) >> 8) + 128;
		}
	}
}



void VideoEncoder::log_packet(const AVFormatContext *fmt_ctx, const AVPacket *pkt)
{
	AVRational *time_base = &fmt_ctx->streams[pkt->stream_index]->time_base;

	char buf1[AV_TS_MAX_STRING_SIZE] = { 0 };
	av_ts_make_string(buf1, pkt->pts);
	char buf2[AV_TS_MAX_STRING_SIZE] = { 0 };
	av_ts_make_string(buf2, pkt->dts);
	char buf3[AV_TS_MAX_STRING_SIZE] = { 0 };
	av_ts_make_string(buf3, pkt->duration);
	char buf4[AV_TS_MAX_STRING_SIZE] = { 0 };
	av_ts_make_time_string(buf4, pkt->pts, time_base);
	char buf5[AV_TS_MAX_STRING_SIZE] = { 0 };
	av_ts_make_time_string(buf5, pkt->dts, time_base);
	char buf6[AV_TS_MAX_STRING_SIZE] = { 0 };
	av_ts_make_time_string(buf6, pkt->duration, time_base);

	printf("pts:%s pts_time:%s dts:%s dts_time:%s duration:%s duration_time:%s stream_index:%d\n",
		buf1, buf4,
		buf2, buf5,
		buf3, buf6,
		pkt->stream_index);
}

int VideoEncoder::write_frame(AVFormatContext *fmt_ctx, const AVRational *time_base, AVStream *st, AVPacket *pkt)
{
	/* rescale output packet timestamp values from codec to stream timebase */
	av_packet_rescale_ts(pkt, *time_base, st->time_base);
	pkt->stream_index = st->index;

	/* Write the compressed frame to the media file. */
	log_packet(fmt_ctx, pkt);
	return av_interleaved_write_frame(fmt_ctx, pkt);
}

/* Add an output stream. */
void VideoEncoder::add_stream(OutputStream *ost, AVFormatContext *oc, AVCodec **codec, enum AVCodecID codec_id, int width, int height)
{
	AVCodecContext *c;
	int i;

	/* find the encoder */
	*codec = avcodec_find_encoder(codec_id);
	if (!(*codec)) {
		fprintf(stderr, "Could not find encoder for '%s'\n",
			avcodec_get_name(codec_id));
		exit(1);
	}

	ost->st = avformat_new_stream(oc, *codec);
	if (!ost->st) {
		fprintf(stderr, "Could not allocate stream\n");
		exit(1);
	}
	ost->st->id = oc->nb_streams - 1;
	c = ost->st->codec;


	c->codec_id = codec_id;

	c->bit_rate = 300000000;
	/* Resolution must be a multiple of two. */
	c->width = width;
	c->height = height;
	/* timebase: This is the fundamental unit of time (in seconds) in terms
	* of which frame timestamps are represented. For fixed-fps content,
	* timebase should be 1/framerate and timestamp increments should be
	* identical to 1. */
	//ost->st->time_base = (AVRational) { 1, STREAM_FRAME_RATE };
	ost->st->time_base.num = 1;
	ost->st->time_base.den = 60;// (AVRational) { 1, c->sample_rate };
	c->time_base = ost->st->time_base;

	c->gop_size = 1; /* emit one intra frame every twelve frames at most */
	c->pix_fmt = AV_PIX_FMT_YUV420P;
	if (c->codec_id == AV_CODEC_ID_MPEG2VIDEO) {
		/* just for testing, we also add B frames */
		c->max_b_frames = 2;
	}

	if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
		/* Needed to avoid using macroblocks in which some coeffs overflow.
		* This does not happen with normal video, it just happens here as
		* the motion of the chroma plane does not match the luma plane. */
		c->mb_decision = 2;
	}

	/* Some formats want stream headers to be separate. */
	if (oc->oformat->flags & AVFMT_GLOBALHEADER)
		c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
}

/**************************************************************/
/* video output */

AVFrame *VideoEncoder::alloc_picture(enum AVPixelFormat pix_fmt, int width, int height)
{
	AVFrame *picture;
	int ret;

	picture = av_frame_alloc();
	if (!picture)
		return NULL;

	picture->format = pix_fmt;
	picture->width = width;
	picture->height = height;

	/* allocate the buffers for the frame data */
	ret = av_frame_get_buffer(picture, 32);
	if (ret < 0) {
		fprintf(stderr, "Could not allocate frame data.\n");
		exit(1);
	}

	return picture;
}

void VideoEncoder::open_video(AVFormatContext *oc, AVCodec *codec, OutputStream *ost, AVDictionary *opt_arg)
{
	int ret;
	AVCodecContext *c = ost->st->codec;
	AVDictionary *opt = NULL;

	av_dict_copy(&opt, opt_arg, 0);

	/* open the codec */
	ret = avcodec_open2(c, codec, &opt);
	av_dict_free(&opt);
	if (ret < 0) {
		char buf1[AV_TS_MAX_STRING_SIZE] = { 0 };
		fprintf(stderr, "Could not open video codec: %s\n", av_make_error_string(buf1, AV_ERROR_MAX_STRING_SIZE, ret));
		exit(1);
	}

	/* allocate and init a re-usable frame */
	ost->frame = alloc_picture(c->pix_fmt, c->width, c->height);
	if (!ost->frame) {
		fprintf(stderr, "Could not allocate video frame\n");
		exit(1);
	}

	/* If the output format is not YUV420P, then a temporary YUV420P
	* picture is needed too. It is then converted to the required
	* output format. */
	ost->tmp_frame = NULL;
	if (c->pix_fmt != AV_PIX_FMT_YUV420P) {
		ost->tmp_frame = alloc_picture(AV_PIX_FMT_YUV420P, c->width, c->height);
		if (!ost->tmp_frame) {
			fprintf(stderr, "Could not allocate temporary picture\n");
			exit(1);
		}
	}
}

/* Prepare a dummy image. */
GLubyte* frame_data = NULL;
void VideoEncoder::fill_yuv_image(AVFrame *pict, int frame_index, int width, int height)
{
	int x, y, i, ret;

	/* when we pass a frame to the encoder, it may keep a reference to it
	* internally;
	* make sure we do not overwrite it here
	*/
	ret = av_frame_make_writable(pict);
	if (ret < 0)
		exit(1);

	i = frame_index;

	if (!frame_data)
	{
		unsigned int size = width * height;
		frame_data = new GLubyte[size * 3];
		for (unsigned int i = 0; i < size; ++i)
		{
			frame_data[i * 3] = 255;
			frame_data[i * 3 + 1] = 128;
			frame_data[i * 3 + 2] = 0;
		}
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, 0);
	glReadBuffer(GL_FRONT);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, frame_data);

	RGB2Yuv420p(pict, frame_data, width, height);


	// Y 
	/*for (y = 0; y < height; y++)
	for (x = 0; x < width; x++)
	pict->data[0][y * pict->linesize[0] + x] = x + y + i * 3;

	// Cb and Cr
	for (y = 0; y < height / 2; y++) {
	for (x = 0; x < width / 2; x++) {
	pict->data[1][y * pict->linesize[1] + x] = 128 + y + i * 2;
	pict->data[2][y * pict->linesize[2] + x] = 64 + x + i * 5;
	}
	}*/
}

AVFrame *VideoEncoder::get_video_frame(OutputStream *ost)
{
	AVCodecContext *c = ost->st->codec;

	/* check if we want to generate more frames */
	AVRational av;
	av.num = 1;
	av.den = 1;

	if (av_compare_ts(ost->next_pts, ost->st->codec->time_base,
		STREAM_DURATION, av) >= 0)
		return NULL;

	if (c->pix_fmt != AV_PIX_FMT_YUV420P) {
		/* as we only generate a YUV420P picture, we must convert it
		* to the codec pixel format if needed */
		if (!ost->sws_ctx) {
			ost->sws_ctx = sws_getContext(c->width, c->height,
				AV_PIX_FMT_YUV420P,
				c->width, c->height,
				c->pix_fmt,
				SCALE_FLAGS, NULL, NULL, NULL);
			if (!ost->sws_ctx) {
				fprintf(stderr,
					"Could not initialize the conversion context\n");
				exit(1);
			}
		}
		fill_yuv_image(ost->tmp_frame, ost->next_pts, c->width, c->height);
		sws_scale(ost->sws_ctx,
			(const uint8_t * const *)ost->tmp_frame->data, ost->tmp_frame->linesize,
			0, c->height, ost->frame->data, ost->frame->linesize);
	}
	else {
		fill_yuv_image(ost->frame, ost->next_pts, c->width, c->height);
	}

	ost->frame->pts = ost->next_pts++;

	return ost->frame;
}

/*
* encode one video frame and send it to the muxer
* return 1 when encoding is finished, 0 otherwise
*/
int VideoEncoder::write_video_frame(AVFormatContext *oc, OutputStream *ost)
{
	int ret;
	AVCodecContext *c;
	AVFrame *frame;
	int got_packet = 0;
	AVPacket pkt = { 0 };

	c = ost->st->codec;

	frame = get_video_frame(ost);

	av_init_packet(&pkt);

	/* encode the image */
	ret = avcodec_encode_video2(c, &pkt, frame, &got_packet);
	if (ret < 0) {
		char buf1[AV_TS_MAX_STRING_SIZE] = { 0 };
		fprintf(stderr, "Error encoding video frame: %s\n", av_make_error_string(buf1, AV_ERROR_MAX_STRING_SIZE, ret));
		system("pause");
		exit(1);
	}

	if (got_packet) {
		ret = write_frame(oc, &c->time_base, ost->st, &pkt);
	}
	else {
		ret = 0;
	}

	if (ret < 0) {
		char buf1[AV_TS_MAX_STRING_SIZE] = { 0 };
		fprintf(stderr, "Error while writing video frame: %s\n", av_make_error_string(buf1, AV_ERROR_MAX_STRING_SIZE, ret));
		system("pause");
		exit(1);
	}

	return (frame || got_packet) ? 0 : 1;
}

void VideoEncoder::close_stream(AVFormatContext *oc, OutputStream *ost)
{
	avcodec_close(ost->st->codec);
	av_frame_free(&ost->frame);
	av_frame_free(&ost->tmp_frame);
	sws_freeContext(ost->sws_ctx);
	swr_free(&ost->swr_ctx);
}