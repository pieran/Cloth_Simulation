#version 150 core

uniform mat4 projViewMatrix;
uniform float z_offset = 10.0f;

in  vec4 position;
in  vec4 colour;

out Vertex {
	vec4 colour;
	vec4 pos;	
} OUT;




void main(void)	{
	vec4 vp = projViewMatrix * vec4(position.xyz, 1.0f);
	
	float linear_z = 0.02 / (1000.01 - vp.z * 999.99);
	linear_z -= z_offset;
	vp.z = ((0.02 / linear_z) - 1000.01) / 999.99;
	gl_Position	  = vp;
	
	OUT.pos		  = vec4(vp.xyz, position.w);	
	OUT.colour    = colour;
}