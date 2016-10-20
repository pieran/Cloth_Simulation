
#pragma once
#include <nclgl\OGLRenderer.h>

class Surface {
public:
	Surface(GLuint width, GLuint height);
	~Surface();

	GLuint GetWidth() { return width; }
	GLuint GetHeight() { return height; }

	void BindSurfaceForDrawing();


	GLuint GetFrameBuffer() { return fbo_surface; }

	int GetNumColorTargets() { return (color_texture1 != 0) ? 2 : 1; }

	void SetColorTexture0(GLuint tex, GLenum type = GL_TEXTURE_2D) { if (color_texture0 == 0) { n_draw_buffers++; }color_texture0 = tex; color_texture0_type = type; }
	void SetColorTexture1(GLuint tex, GLenum type = GL_TEXTURE_2D) { if (color_texture1 == 0) { n_draw_buffers++; } color_texture1 = tex; color_texture1_type = type; }
	void SetDepthStencilTexture0(GLuint tex, bool useDepthComponent = true, bool useStencilComponent = true, GLenum type = GL_TEXTURE_2D) {
		depth_texture0 = tex;
		depth_texture0_type = type;
		use_depth_component = useDepthComponent;
		use_stencil_component = useStencilComponent;
	}
	void SetDepthTexture0(GLuint tex, GLenum type = GL_TEXTURE_2D) {
		depth_texture0 = tex;
		depth_texture0_type = type;
		use_depth_component = true;
	}
	void SetStencilTexture0(GLuint tex, GLenum type = GL_TEXTURE_2D) {
		depth_texture0 = tex;
		depth_texture0_type = type;
		use_stencil_component = true;
	}

	void SetDepthStencilEnabled(bool useDepthComponent = true, bool useStencilComponent = true) {
		use_depth_component = useDepthComponent;
		use_stencil_component = useStencilComponent;
	}

	GLuint GetColorTexture(int index) { return (index == 0) ? color_texture0 : color_texture1; }
	GLuint GetColorTexture0() { return color_texture0; }
	GLuint GetColorTexture1() { return color_texture1; }
	GLuint GetDepthTexture0() { return depth_texture0; }
	GLuint GetStencilTexture0() { return depth_texture0; }

	void CreateColorTexture0(GLint internalFormat = GL_RGB8, bool linear = false, bool repeat = false);
	void CreateColorTexture1(GLint internalFormat = GL_RGB8, bool linear = false, bool repeat = false);
	void CreateDepthTexture0(GLint internalFormat = GL_DEPTH24_STENCIL8, bool linear = false, bool repeat = false);

	bool GenerateFBO();

	void Resize(int width, int height);
	void SetDrawBuffers(unsigned int num_buffers, const GLenum* buffersunsigned);
protected:
	void BuildTextureComponent(GLuint* texture, GLint fbo_component, GLint internalFormat, bool linear, bool repeat);

protected:
	GLuint width, height;

	GLuint color_texture0;
	GLuint color_texture1;
	GLuint depth_texture0;
	bool use_depth_component, use_stencil_component;

	GLenum color_texture0_type;
	GLenum color_texture1_type;
	GLenum depth_texture0_type;

	GLuint n_draw_buffers;
	GLenum draw_buffers[2];


	GLuint fbo_surface;
};