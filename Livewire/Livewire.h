#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_Livewire.h"

class MyInteraction;
class MyGLCanvas;
class PathUtils;

class Livewire : public QMainWindow
{
	Q_OBJECT

public:
	Livewire(QWidget *parent = Q_NULLPTR);

private:
	Ui::LivewireClass ui;

public slots:

	void OpenMeshFile();
	void SaveCurveFile();

	void Geodesy();
	void AnisGeodesyMin();
	void AnisGeodesyMax();

public:

	MyInteraction* m_interaction;
	PathUtils* m_pathUtils;

private:

	void ResetAll();
	void ConstructMesh(const char* fileName, const char* fileType);

	void keyPressEvent(QKeyEvent *e);
	void keyReleaseEvent(QKeyEvent *e);
	void wheelEvent(QWheelEvent *e);
	void resizeEvent(QResizeEvent* event);
};
