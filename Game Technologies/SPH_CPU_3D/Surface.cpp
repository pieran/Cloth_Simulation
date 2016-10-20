
#include "Surface.h"


Surface::Surface(GLuint width, GLuint height) {
	this->width = width;
	this->height = height;

	color_texture0_type = GL_TEXTURE_2D;
	color_texture1_type = GL_TEXTURE_2D;
	depth_texture0_type = GL_TEXTURE_2D;

	fbo_surface = NULL;
	color_texture0 = NULL;
	color_texture1 = NULL;
	depth_texture0 = NULL;

	use_depth_component = true;
	use_stencil_component = true;

	n_draw_buffers = 0;
	draw_buffers[0] = GL_COLOR_ATTACHMENT0;
	draw_buffers[1] = GL_COLOR_ATTACHMENT1;
}

Surface::~Surface() {
	glDeleteTextures(1, &color_texture0);
	glDeleteTextures(1, &color_texture1);
	glDeleteTextures(1, &depth_texture0);

	glDeleteFramebuffers(1, &fbo_surface);
}

void Surface::BindSurfaceForDrawing() {
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_surface);
	glViewport(0, 0, width, height);
}



void Surface::BuildTextureComponent(GLuint* texture, GLint fbo_component, GLint internalFormat, bool linear, bool repeat) {
	if (*texture == NULL)
		glGenTextures(1, texture);

	glBindTexture(GL_TEXTURE_2D, *texture);
	glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, width, height, 0,
		(fbo_component == GL_DEPTH_ATTACHMENT) ? GL_DEPTH_COMPONENT : GL_RGB,
		(fbo_component == GL_DEPTH_ATTACHMENT) ? GL_FLOAT : GL_UNSIGNED_BYTE, NULL);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, linear ? GL_LINEAR : GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear ? GL_LINEAR : GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, repeat ? GL_REPEAT : GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, repeat ? GL_REPEAT : GL_CLAMP_TO_EDGE);

	if (fbo_component == GL_DEPTH_ATTACHMENT) {

		//Has to be GL_NONE in order to correctly sample from glsl shader, 
		//	enable again to allow for hardware optimisation (but without the chance to correctly sample manually)
		//glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);


		//Handle depth as a normal texture/sampler2D
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	}

	glBindTexture(GL_TEXTURE_2D, 0);
}


void Surface::CreateColorTexture0(GLint internalFormat, bool linear, bool repeat) {
	if (color_texture0 == NULL) n_draw_buffers++;

	BuildTextureComponent(&color_texture0, GL_COLOR_ATTACHMENT0, internalFormat, linear, repeat);
}
void Surface::CreateColorTexture1(GLint internalFormat, bool linear, bool repeat) {
	if (color_texture1 == NULL) n_draw_buffers++;

	BuildTextureComponent(&color_texture1, GL_COLOR_ATTACHMENT1, internalFormat, linear, repeat);
}
void Surface::CreateDepthTexture0(GLint internalFormat, bool linear, bool repeat) {
	BuildTextureComponent(&depth_texture0, GL_DEPTH_ATTACHMENT, internalFormat, linear, repeat);
	use_stencil_component = (internalFormat == GL_DEPTH24_STENCIL8);
	use_depth_component = true;
}

bool Surface::GenerateFBO() {
	if (n_draw_buffers > 2) n_draw_buffers = 2;

	if (fbo_surface == NULL) {
		glGenFramebuffers(1, &fbo_surface);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, fbo_surface);

	if (depth_texture0 != NULL) {
		if (use_depth_component) {
			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_texture0_type, depth_texture0, 0);
		}
		if (use_stencil_component) {
			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, depth_texture0_type, depth_texture0, 0);
		}
	}

	if (color_texture0 != NULL) {
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, color_texture0_type, color_texture0, 0);
	}

	if (color_texture1 != NULL) {
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, color_texture1_type, color_texture1, 0);
	}

	glDrawBuffers(n_draw_buffers, draw_buffers);


	GLenum status;
	status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE) {
		printf("FRAMEBUFFER ERROR: %d\n", status);
		return false;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	return true;
}

void Surface::SetDrawBuffers(unsigned int num_buffers, const GLenum* buffers)
{
	if (num_buffers <= 2)
	{
		n_draw_buffers = num_buffers;
		memcpy(draw_buffers, buffers, n_draw_buffers * sizeof(GLenum));

		glBindFramebuffer(GL_FRAMEBUFFER, fbo_surface);
		glDrawBuffers(n_draw_buffers, draw_buffers);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
}