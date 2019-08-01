#include "Livewire.h"

#include <QtGui/QtGui>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <iostream>
#include <QtWidgets/QInputDialog>

#include "core/Interaction.h"
#include "core/PathUtils.h"
#include "MyGLCanvas.h"

Livewire::Livewire(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	//before everything get started, initialize the mesh, interaction and render tool.
	m_pathUtils = new PathUtils();
	*(ui.openGLWidget->getPathUtils()) = m_pathUtils;

	m_interaction = new MyInteraction();
	*(ui.openGLWidget->getInteractionUtils()) = m_interaction;
	*(m_pathUtils->getInteractionUtils()) = m_interaction;

	//ui.openGLWidget->initializeGL();
	ui.openGLWidget->pos_y = this->frameGeometry().height() - ui.openGLWidget->frameGeometry().height() - ui.menuBar->frameGeometry().height();

	ui.Geodesy_checkBox->setVisible(false);
	ui.AnisGeoMin_checkBox->setVisible(false);
	ui.AnisGeoMax_checkBox->setVisible(false);
}

void Livewire::ResetAll()
{
	delete m_pathUtils;
	delete m_interaction;

	m_pathUtils = new PathUtils();
	*(ui.openGLWidget->getPathUtils()) = m_pathUtils;

	m_interaction = new MyInteraction();
	*(m_pathUtils->getInteractionUtils()) = m_interaction;
	*(ui.openGLWidget->getInteractionUtils()) = m_interaction;

	ui.openGLWidget->initializeGL();
	ui.openGLWidget->setMouseTracking(true);
}

//first step: load a mesh, with extension of off,obj,ply...
void Livewire::OpenMeshFile()
{
	QString file_name = QFileDialog::getOpenFileName(this,
		tr("Open File"),
		"",
		"*.off;*.ply;*.obj;;*.off;;*.ply;;*.obj",
		0);

	if (!file_name.isNull())
	{
		char tFileName[400];
		strcpy(tFileName, file_name.toStdString().data());
		char FileType[400];
		for (int i = strlen(tFileName) - 1; i >= 0; i--) {
			if (tFileName[i] == '.') {
				i++;
				unsigned tsize = strlen(tFileName) - i;
				for (unsigned j = 0; j < tsize; j++) {
					FileType[j] = tFileName[i + j];
				}
				FileType[tsize] = '\0';
				break;
			}
		}

		ConstructMesh(file_name.toStdString().data(), FileType);
		ui.openGLWidget->paintGL();

		//OverSegmentation();
		//ui.openGLWidget->paintGL();

		ui.Geodesy_checkBox->setVisible(true);
		ui.AnisGeoMin_checkBox->setVisible(true);
		ui.AnisGeoMax_checkBox->setVisible(true);
	}
}
void Livewire::ConstructMesh(const char* fileName, const char* fileType)
{
	ResetAll();

	m_pathUtils->m_filename = fileName;

	if (strcmp(fileType, "ply") == 0 || strcmp(fileType, "off") == 0 || strcmp(fileType, "obj") == 0)
	{
		m_pathUtils->readTrimesh(fileName);
	}
	else
		return;

	*(m_interaction->getMyMesh()) = *(m_pathUtils->getMyMesh());
	m_interaction->initSelectionColor();

	//m_pathUtils->globalParameter();

	srand(time(NULL));
}

void Livewire::SaveCurveFile()
{
	m_pathUtils->SavePath();

	//QMessageBox message(QMessageBox::NoIcon, "Save File", QString("The files are save to the same directory. There are three files, including a per vertex labeling file(.label), two segmentation boundary files(.curve), with or without smoothing."));
	//message.exec();
}

void Livewire::Geodesy()
{
	m_pathUtils->metricActived(EUCLIDEAN, ui.Geodesy_checkBox->isChecked());
	ui.openGLWidget->paintGL();
}
void Livewire::AnisGeodesyMin()
{
	m_pathUtils->metricActived(OUR_MIN, ui.AnisGeoMin_checkBox->isChecked());
	ui.openGLWidget->paintGL();
}
void Livewire::AnisGeodesyMax()
{
	m_pathUtils->metricActived(OUR_MAX, ui.AnisGeoMax_checkBox->isChecked());
	ui.openGLWidget->paintGL();
}


void Livewire::keyPressEvent(QKeyEvent *e)
{

	ui.openGLWidget->keyPressEvent(e);
}
void Livewire::keyReleaseEvent(QKeyEvent *e)
{

	ui.openGLWidget->keyReleaseEvent(e);
}

void Livewire::wheelEvent(QWheelEvent *e)
{
	//ui.openGLWidget->wheelEvent(e);

}

void Livewire::resizeEvent(QResizeEvent* event)
{
	QMainWindow::resizeEvent(event);
	// Your code here
	//int width = this->frameGeometry().width();
	//int height = this->frameGeometry().height();
	int width = ui.centralWidget->frameGeometry().width();
	int height = ui.centralWidget->frameGeometry().height();
	ui.openGLWidget->setFixedSize(QSize(width, height));

	ui.openGLWidget->resizeGL(width, height);
	ui.openGLWidget->paintGL();
}