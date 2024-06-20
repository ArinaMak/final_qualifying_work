#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "materials.h"
#include "dialogchangeparam.h"
// #include "formchangeparam.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    QColor color_graph=Qt::red;
    double rho;
    double moduleY;
    Materials mat123;
    int number_of_subtask=0;


private slots:
    bool isInputCorrect(QString &str);
    double stringToNumber(QString &str, QString &str_if_not_correct, QString str_name);
    void plotGraphByPoints(QVector<double> &vec_x, QVector<double> &vec_y, bool isDots);

    QString printingInformationAboutGraphs();
    void fillingInDataFromString(QVector<double> &vec_x, QString &str);
    void OpenTxtFile();
    void on_actionOpen_triggered();
    void on_pushButton_PlotGraph_clicked();
    void on_pushButton_ClearGraph_clicked();

    void on_pushButton_Go_clicked();

    void on_comboBox_Material_activated(int index);
    void on_pushButton_ColorChange_clicked();

    void on_pushButton_Move_clicked(bool checked);
    void on_pushButton_ClearGraph_pressed();
    void on_pushButton_SaveGraph_clicked();
    void on_pushButton_SaveImage_clicked();
    void on_pushButton_OpenFile_clicked();
    void on_pushButton_Zoom_clicked(bool checked);

    void on_pushButton_AddDot_clicked();

    void on_pushButton_ExampleTask_clicked();

    void on_pushButton_SaveMaterial_clicked();



private:
    Ui::MainWindow *ui;

    DialogChangeParam *dialog;
    // FormChangeParam *form;

public slots:
    void slotForm(int number);
};
#endif // MAINWINDOW_H
