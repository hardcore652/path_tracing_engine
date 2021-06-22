#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>

using namespace sf;
using namespace std;

int windowWidth = 1600;
int windowHeight = 900;

float skyBrightness = 1.0f;
float mouseSensitivity = 2.0;
float cameraAperture = 0.05;
float focalDistance = 5.0f;
bool autoFocalDistance = true;
int framerate = 0;

void clearConsole() {
	cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
}
void printDebug() {
	cout << "FPS: " << framerate << endl;
	cout << "Depth of field:\n	Focal Distance: ";
	if (autoFocalDistance) cout << "Auto";
	else cout << focalDistance;
	cout << "\n	Aperture: " << cameraAperture << endl;
	cout << "Sky brightness: " << skyBrightness << endl;
}

int main() {
	int mouseX = 0;
	int mouseY = 0;
	bool firstPersonMode = true;
	int frame = 1;
	float dt = 0.0;
	float previousTime = 0.0;
	float time = 0.0;
	int framesCounter = 0;
	float startFramesCounterTime = 0;
	Vector3f cameraPos = Vector3f(-5.0, 0.0, 0.0);
	float cameraSpeed = 0.2;
	bool movement[6] = { false, false, false, false, false, false };
	bool altHotkeys = false;
	float debugDuration = 0.0;

	Clock clock; // Clock object
	system("mode 40, 8");
	RenderWindow window(VideoMode(windowWidth, windowHeight), "hcPathTracing"); // Create the window
	window.setMouseCursorVisible(false);

	const String shaderPath = "shader.frag";
	const String skyboxPath = "data\\skybox.jpg";

	Shader shader;
	shader.loadFromFile(shaderPath, Shader::Fragment);
	Texture skybox;
	skybox.loadFromFile(skyboxPath);
	shader.setUniform("u_resolution", Vector2f(windowWidth, windowHeight));
	shader.setUniform("u_skybox", skybox);
	shader.setUniform("u_skyBrightness", skyBrightness);
	shader.setUniform("u_cameraAperture", cameraAperture);

	RenderTexture firstBuff; // we will use doublebuff
	RenderTexture secondBuff;
	firstBuff.create(windowWidth, windowHeight);
	secondBuff.create(windowWidth, windowHeight);
	Sprite firstBuffSprite = Sprite(firstBuff.getTexture());
	Sprite secondBuffSprite = Sprite(secondBuff.getTexture());
	Sprite firstBuffFlippedSprite = Sprite(firstBuff.getTexture()); // Flipped vertical rendering order of pixels
	Sprite secondBuffFlippedSprite = Sprite(secondBuff.getTexture());
	firstBuffFlippedSprite.setScale(1, -1); // Flip and move sprites
	firstBuffFlippedSprite.setPosition(0, windowHeight);
	secondBuffFlippedSprite.setScale(1, -1);
	secondBuffFlippedSprite.setPosition(0, windowHeight);

	random_device rd;
	mt19937 e2(rd());
	uniform_real_distribution<> dist(0.0f, 1.0f);

	
	while (window.isOpen()) {

		Event event;
		while (window.pollEvent(event)) {
			if (event.type == Event::Closed) window.close();
			else if (event.type == Event::MouseMoved)
			{
				if (firstPersonMode)
				{
					int dx = event.mouseMove.x - windowWidth / 2;
					int dy = event.mouseMove.y - windowHeight / 2;
					mouseX += dx;
					mouseY += dy;
					Mouse::setPosition(Vector2i(windowWidth / 2, windowHeight / 2), window);
					if (dx != 0 || dy != 0) frame = 1;
				}
			}
			else if (event.type == Event::MouseButtonPressed) {
				window.setMouseCursorVisible(false);
				firstPersonMode = true;
			}
			else if (event.type == Event::MouseWheelScrolled) {
				if (event.mouseWheelScroll.wheel == Mouse::VerticalWheel) {
					if (altHotkeys) { cameraAperture += event.mouseWheelScroll.delta / 20.0f; if (cameraAperture < 0.0) cameraAperture = 0.0f; shader.setUniform("u_cameraAperture", cameraAperture); }
					else if (autoFocalDistance) {
						skyBrightness += event.mouseWheelScroll.delta / 5.0f;
						shader.setUniform("u_skyBrightness", skyBrightness);
					}
					else focalDistance += event.mouseWheelScroll.delta / 2.5f;
					frame = 1;
				}
			}
			else if (event.type == Event::KeyPressed)
			{
				if (event.key.code == Keyboard::Escape)
				{
					window.setMouseCursorVisible(true);
					firstPersonMode = false;
				}
				else if (event.key.code == Keyboard::R) {
					shader.loadFromFile(shaderPath, Shader::Fragment);
					skybox.loadFromFile(skyboxPath);
					shader.setUniform("u_resolution", Vector2f(windowWidth, windowHeight));
					shader.setUniform("u_skybox", skybox);
					shader.setUniform("u_skyBrightness", skyBrightness);
					shader.setUniform("u_cameraAperture", cameraAperture);
					frame = 1;
				}
				else if (event.key.code == Keyboard::LAlt) altHotkeys = true;
				else if (event.key.code == Keyboard::Q) {
					if (altHotkeys) window.close();
				}
				else if (altHotkeys && event.key.code == Keyboard::D) {
					autoFocalDistance = not autoFocalDistance;
					frame = 1;
				}
				else if (not altHotkeys && (event.key.code == Keyboard::W || event.key.code == Keyboard::A || event.key.code == Keyboard::S || event.key.code == Keyboard::D || event.key.code == Keyboard::LShift || event.key.code == Keyboard::Space)) {
					window.setMouseCursorVisible(false);
					firstPersonMode = true;
					if (event.key.code == Keyboard::W) movement[0] = true;
					else if (event.key.code == Keyboard::A) movement[2] = true;
					else if (event.key.code == Keyboard::S) movement[1] = true;
					else if (event.key.code == Keyboard::D) movement[3] = true;
					else if (event.key.code == Keyboard::LShift) movement[4] = true;
					else if (event.key.code == Keyboard::Space) movement[5] = true;
				}
			}
			else if (event.type == Event::KeyReleased)
			{
				if (event.key.code == Keyboard::W) movement[0] = false;
				else if (event.key.code == Keyboard::A) movement[2] = false;
				else if (event.key.code == Keyboard::S) movement[1] = false;
				else if (event.key.code == Keyboard::D) movement[3] = false;
				else if (event.key.code == Keyboard::LShift) movement[4] = false;
				else if (event.key.code == Keyboard::Space) movement[5] = false;
				else if (event.key.code == Keyboard::LAlt) altHotkeys = false;
			}
		}

		time = clock.getElapsedTime().asSeconds();
		dt = time - previousTime;
		dt *= 60;
		float s = cameraSpeed * dt;
		float mx = (float(mouseX) / windowWidth) * mouseSensitivity;
		float my = (float(mouseY) / windowHeight) * mouseSensitivity;
		if (my > 1.57) {
			my = 1.57;
			mouseY = my / mouseSensitivity * windowHeight;
		}
		else if (my < -1.57) {
			my = -1.57;
			mouseY = my / mouseSensitivity * windowHeight;
		}
		if (movement[0]) {
			cameraPos.x += s * cos(mx);
			cameraPos.y += s * sin(mx);
		}
		if (movement[1]) {
			cameraPos.x -= s * cos(mx);
			cameraPos.y -= s * sin(mx);
		}
		if (movement[2]) {
			cameraPos.x += s * sin(mx);
			cameraPos.y -= s * cos(mx);
		}
		if (movement[3]) {
			cameraPos.x -= s * sin(mx);
			cameraPos.y += s * cos(mx);
		}
		if (movement[4]) cameraPos.z += s;
		if (movement[5]) cameraPos.z -= s;
		for (int i = 0; i < 6; i++) {
			if (movement[i]) {
				frame = 1; break;
			}
		}

		shader.setUniform("u_time", time);
		shader.setUniform("u_mouse", Vector2f(mx, my));
		shader.setUniform("u_cameraPos", cameraPos);
		shader.setUniform("u_samplePart", 1.0f / frame);
		shader.setUniform("u_seed", Vector2f((float)dist(e2), (float)dist(e2)) * 999.0f);
		shader.setUniform("u_seed2", Vector2f((float)dist(e2), (float)dist(e2)) * 999.0f);
		if (autoFocalDistance) shader.setUniform("u_focalDist", -1.0f);
		else shader.setUniform("u_focalDist", focalDistance);
		if (time - debugDuration >= 0.3) {
			clearConsole();
			printDebug();
			debugDuration = time;
		}
		if (time - startFramesCounterTime >= 1.0) {
			framerate = framesCounter;
			framesCounter = 0;
			startFramesCounterTime = time;
		}

		if (frame % 2 == 1) {
			shader.setUniform("u_tex", firstBuff.getTexture());
			secondBuff.draw(firstBuffFlippedSprite, &shader);
			window.draw(secondBuffSprite);
		}
		else {
			shader.setUniform("u_tex", secondBuff.getTexture());
			firstBuff.draw(secondBuffFlippedSprite, &shader);
			window.draw(firstBuffSprite);
		}

		window.display();
		previousTime = time;
		frame++;
		framesCounter++;
	}
	return 0;
}