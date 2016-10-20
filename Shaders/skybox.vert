#version 150 core

in vec3 position;
uniform mat4 projViewMatrix;

out Vertex	{
	vec3 texCoords;
} OUT;


void main(void)	{
	OUT.texCoords = position;
	gl_Position = projViewMatrix * vec4(position.xyz,1.0f);
}