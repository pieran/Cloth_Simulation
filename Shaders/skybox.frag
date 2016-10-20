#version 150 core

uniform samplerCube texCubemap;

in Vertex	{
	vec3 texCoords;
} IN;

out vec4 gl_FragColor;


void main(void)	{
	//vec3 tcoords = invViewMatrix * (vec3(IN.texCoords.x, IN.texCoords.y, 0.0f) * 2.0f - 1.0f);
	gl_FragColor 	= texture(texCubemap, IN.texCoords);
}