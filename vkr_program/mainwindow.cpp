#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <QVector>
#include "material.h"
#include "materials.h"
#include "task.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    // form=new FormChangeParam;
    dialog = new DialogChangeParam;

    //connect(dialog, &DialogChangeParam::signalForm, this, &MainWindow::slotForm);
    //connect(form, &FormChangeParam::signalForm, this, &MainWindow::slotForm);
    connect(dialog, &DialogChangeParam::signalForm, this, &MainWindow::slotForm);


    QString s0="background: #ff0000;\nborder-radius: 5px;\nborder: 3px solid white;";

    ui->pushButton_ColorChange->setStyleSheet(s0);
    ui->pushButton_ColorChange->update();
    ui->widget_Graph->xAxis->setLabel("x");
    ui->widget_Graph->yAxis->setLabel("F");

    QString s1="background: #eaeaea;\nborder-radius: 5px;\nborder: 3px solid #eaeaea;";
    ui->pushButton_Move->setStyleSheet(s1);
    ui->pushButton_Move->update();

    ui->pushButton_Zoom->setStyleSheet(s1);
    ui->pushButton_Zoom->update();

    ui->pushButton_ClearGraph->setStyleSheet(s1);
    ui->pushButton_ClearGraph->update();

    ui->widget_Graph->addGraph();   // Добавление нулевого графа для пользовательских точек

    ui->lineEdit_F_1->setPlaceholderText(QString::fromUtf8("м\u{00b2}"));
    ui->lineEdit_F_2->setPlaceholderText(QString::fromUtf8("м\u{00b2}"));
    // ui->pushButton_Go->setStyleSheet("background-color: rgba(255, 255, 255, 0);");
    // ui->pushButton_Save->setStyleSheet("background-color: rgba(255, 255, 255, 0);");


    ui->comboBox_Material->addItem("Ввести другой материал");
    for(int i=0; i<mat123.material_list.length(); i++)
    {
        ui->comboBox_Material->addItem(mat123.material_list.value(i)->aboutMaterial());
    }
    // ui->comboBox_Material->view()->setMinimumWidth(7*mat123.max_length);
    ui->comboBox_Material->view()->setMinimumWidth(14*mat123.max_length);
    ui->comboBox_Material->setCurrentIndex(-1);
    ui->comboBox_Material->setToolTip(ui->comboBox_Material->currentText());

    ui->label_MaterialName->setVisible(false);
    ui->label_MaterialDensity->setVisible(false);
    ui->label_MaterialYModule->setVisible(false);

    ui->lineEdit_MaterialName->setVisible(false);
    ui->lineEdit_MaterialDensity->setVisible(false);
    ui->lineEdit_MaterialYModule->setVisible(false);

    ui->pushButton_SaveMaterial->setVisible(false);

    QString s("background: #"
              + QString(color_graph.red() < 16? "0" : "") + QString::number(color_graph.red(),16)
              + QString(color_graph.green() < 16? "0" : "") + QString::number(color_graph.green(),16)
              + QString(color_graph.blue() < 16? "0" : "") + QString::number(color_graph.blue(),16)
              + ";\nborder-radius: 5px;"
              + "\nborder: 3px solid white;");
    ui->pushButton_ColorChange->setStyleSheet(s);
    ui->pushButton_ColorChange->update();

    ui->comboBox_TypeOfTask->view()->setMinimumWidth(300);
}

MainWindow::~MainWindow()
{
    ui->comboBox_Material->clear();
    ui->widget_Graph->clearGraphs();
    // delete form;
    delete dialog;
    delete ui;
}

bool MainWindow::isInputCorrect(QString &str)
{
    int length=str.length();
    char index_e=0;
    char index_dot=0;
    char index_sign=0;
    for(int i=0; i<length; i++)
    {
        if(!str[i].isDigit())
        {
            if(str[i]=='.')
            {
                if(index_dot!=0) return false;
                index_dot=i;
                if(index_e!=0 && index_dot>index_e) return false;
                if(index_dot+1<length && str[index_dot+1].isDigit()) i++;
                else return false;
                continue;
            }
            if(str[i]=='e')
            {
                if(index_e!=0) return false;
                index_e=i;
                if(index_dot!=0 && index_dot>index_e) return false;

                if(index_e+1<length)
                {
                    i++;
                    if(str[i]=="+" || str[i]=="-")
                    {
                        if(i+1<length && str[i+1].isDigit()) i++;
                        else return false;
                        continue;
                    }
                    if(str[i].isDigit()) i++;
                    else return false;
                }
                else return false;
                continue;
            }
            return false;
        }
    }
    return true;
}

double MainWindow::stringToNumber(QString &str, QString &str_if_not_correct, QString str_name)   //
{
    double tmp_number=0.0;
    if(isInputCorrect(str)) //
    {
        tmp_number=str.toDouble();
        if(tmp_number==0 && str_name=="массы") return tmp_number;
        else if(tmp_number>0)  return tmp_number;   // не обязательна проверка, т.к. попадают только положительные числа
    }

    if(str_if_not_correct!="") str_if_not_correct+=", ";
    str_if_not_correct+=str_name;

    return 0.0;
}

void MainWindow::plotGraphByPoints(QVector<double> &vec_x, QVector<double> &vec_y, bool isDots)
{
    int last_graph=ui->widget_Graph->graphCount();
    if(isDots==true)
    {
        ui->widget_Graph->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
        ui->widget_Graph->graph(0)->setLineStyle(QCPGraph::lsNone);
        last_graph=0;
    }
    else ui->widget_Graph->addGraph();

    ui->widget_Graph->graph(last_graph)->addData(vec_x, vec_y);
    ui->widget_Graph->graph(last_graph)->setPen(QPen(color_graph,1));
    ui->widget_Graph->replot();

    qDebug()<<"x_i = "<<vec_x;
    qDebug()<<"F_i = "<<vec_y;
}

QString MainWindow::printingInformationAboutGraphs()
{
    int number_of_graphs=ui->widget_Graph->graphCount();
    int curr_graph=0;

    QString info="";
    info+="Graphs = "+ QString::number(number_of_graphs)+":\n";

    for(curr_graph=0; curr_graph<number_of_graphs; curr_graph++)
    {
        info+="["+QString::number(curr_graph)+"]:\n";
        int number_of_values=ui->widget_Graph->graph(curr_graph)->dataCount();
        info+="x_i="+QString::number(number_of_values)+": ";
        QString info_x_i="";
        QString info_F_i="";
        for(int j=0; j<number_of_values-1; j++)
        {
            double key_tmp=ui->widget_Graph->graph(curr_graph)->data()->at(j)->key;
            double value_tmp=ui->widget_Graph->graph(curr_graph)->data()->at(j)->value;
            info_x_i+=QString::number(key_tmp)+", ";
            info_F_i+=QString::number(value_tmp)+", ";
        }
        if(number_of_values>0)
        {
            info_x_i+=QString::number(ui->widget_Graph->graph(curr_graph)->data()->at(number_of_values-1)->key)+";\n";
            info_F_i+=QString::number(ui->widget_Graph->graph(curr_graph)->data()->at(number_of_values-1)->value)+";\n";
        }
        else
        {
            info_x_i+=";\n";
            info_F_i+=";\n";
        }
        info+=info_x_i+"F_i: "+info_F_i+"color: "+QString::number(ui->widget_Graph->graph(curr_graph)->pen().color().red())
                +", "+QString::number(ui->widget_Graph->graph(curr_graph)->pen().color().green())
                +", "+QString::number(ui->widget_Graph->graph(curr_graph)->pen().color().blue()) +";\n";

    }
    qDebug()<<"Имеющаяся информация:\n\n"<<info;
    return info;
}

void MainWindow::fillingInDataFromString(QVector<double> &vec_x, QString &str)
{

    int tmp = str.indexOf(":");
    int tmp2 = str.indexOf(";");

    QString tmp_line=str.mid(tmp+1, tmp2-tmp-1);
    QStringList list_numbers=tmp_line.split(',');

    for(int j=0; j<list_numbers.length(); j++)
    {
        qDebug()<<"цикл = "<<list_numbers[j];
                                         vec_x.push_back((list_numbers[j]).toDouble());
    }
    qDebug()<<"\n\nvec_x="<<vec_x<<"\n";
    qDebug()<<"vector length = "<<vec_x.length();
}

void MainWindow::OpenTxtFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "E:/Programming/tests/test_save", tr("Text Files (*.txt)"));

    QFile fileIn(fileName);

    if(fileIn.open(QIODevice::ReadOnly))
    {
        QTextStream in(&fileIn);
        int number_of_curr_line=0;
        in.readLine();

        QVector<double> x_i, F_i, color_i;
        bool is_first_graph=true;

        while (!in.atEnd())
        {
            QString line = in.readLine();


            if(line[line.length()-1]==':' || line[line.length()-1]==';') number_of_curr_line++;

            if(number_of_curr_line%2==0 && is_first_graph && line.indexOf("=0:")!=-1)
            {
                is_first_graph=false;
                in.readLine();
                in.readLine();
                number_of_curr_line=0;
            }
            else if(number_of_curr_line%2==0 && number_of_curr_line%4==2) fillingInDataFromString(x_i, line);
            else if(number_of_curr_line%3==0) fillingInDataFromString(F_i, line);
            else if(number_of_curr_line%4==0)
            {
                fillingInDataFromString(color_i, line);
                number_of_curr_line=0;

                color_graph=QColor(color_i[0], color_i[1], color_i[2]);
                plotGraphByPoints(x_i, F_i, is_first_graph);
                if(is_first_graph) is_first_graph=false;

                x_i.clear();
                F_i.clear();
                color_i.clear();
            }
        }
        fileIn.close();
    }

}


void MainWindow::on_actionOpen_triggered()
{
    OpenTxtFile();
}

void MainWindow::on_pushButton_PlotGraph_clicked()
{
    if(color_graph==Qt::red)
    {
        QVector<double> x_i={0,1,2,3,4};
        QVector<double> F_i={1,1,1,1,1};
        plotGraphByPoints(x_i, F_i, false);
    }
    else
    {
        qDebug()<<"hello";
        QVector<double> x_i,F_i;
        for(int i=0; i<300; i++)
        {
            x_i.push_back(i);
            F_i.push_back(i);
        }

        plotGraphByPoints(x_i, F_i, false);
        x_i.clear();
        F_i.clear();
    }
}


void MainWindow::on_pushButton_ClearGraph_clicked()
{
    QString s1="background: #eaeaea;\nborder-radius: 5px;\nborder: 3px solid #eaeaea;";
    ui->pushButton_ClearGraph->setStyleSheet(s1);
    ui->pushButton_ClearGraph->update();

    // qDebug()<<"\nДо очистки кол-во графиков = "<<ui->widget_Graph->graphCount();
    ui->widget_Graph->clearGraphs();
    // qDebug()<<"После очистки кол-во графиков = "<<ui->widget_Graph->graphCount()<<"\n";
    ui->widget_Graph->replot();

    ui->widget_Graph->addGraph();   // Добавление нулевого графика для точек
}


void MainWindow::on_pushButton_Go_clicked()
{
    qDebug()<<"click go!";
    QString l_text=ui->lineEdit_l->text();
    QString F_1_text=ui->lineEdit_F_1->text();
    QString F_2_text=ui->lineEdit_F_2->text();
    QString M_text=ui->lineEdit_M->text();
    QString omega_text=ui->lineEdit_omega->text();

    QString str_if_not_correct="";

    QString info="";

    double l = stringToNumber(l_text, str_if_not_correct, "длины");
    double F_1 = stringToNumber(F_1_text, str_if_not_correct, "минимальной площади поперечного сечения (F_1)");
    double F_2 = stringToNumber(F_2_text, str_if_not_correct, "максимальной площади поперечного сечения (F_2)");
    double M = stringToNumber(M_text, str_if_not_correct, "массы");
    double omega = stringToNumber(omega_text, str_if_not_correct, "частоты");

    qDebug()<<"str_if_not_correct after check = "<<str_if_not_correct;

    if(str_if_not_correct!="")
    {
        str_if_not_correct+=". Строки не должны содержать буквы и отрицательные значения.";
        QString info_str="Некорректный ввод данных задачи!\nПроверьте правильность ввода "+str_if_not_correct;
        QMessageBox::warning(this, "Ввод данных", info_str);
    }
    else
    {
        info+="введено\nl = "+QString::number(l)+"\nF_1 = "+QString::number(F_1)+"\nF_2 = "+QString::number(F_2)+"\nM = "+QString::number(M)+"\nomega = "+QString::number(omega)+"\nrho = "+QString::number(rho)+"\nmpdule Y = "+QString::number(moduleY);

        Task task123(l, F_1, F_2, M, omega, rho, moduleY, 0.00001, 100, 1);

        info+="\nbeta = "+QString("%1").arg(task123.beta, 0, 'f', 8)+"\na = "+QString("%1").arg(task123.a, 0, 'f', 8);

        QVector<double> x_i;
        QVector<double> F_i;

        if(ui->comboBox_TypeOfTask->currentIndex()==0)
        {
            if(task123.isTheOptimalPoint_OnePart(task123.l))
            {
                qDebug()<<"123task123.isTheOptimalPoint_OnePart(task123.l) = "<<task123.isTheOptimalPoint_OnePart(task123.l);
                info+="\nВведенные данные являются оптимальными";
                task123.fillingInVectors_OnePart(x_i, F_i);
                plotGraphByPoints(x_i, F_i, false);
            }
            else
            {
                info+="\nВведенные данные не являются оптимальными";
                dialog->setModal(true);
                dialog->exec();

                if(number_of_subtask==1) task123.new_M_and_l_FromSpecifiedF_OnePart();
                else if(number_of_subtask==2) task123.new_l_and_F2_FromSpecifiedF1_and_M_OnePart();
                else if(number_of_subtask==3) task123.new_l_and_F1_FromSpecifiedF2_and_M_OnePart();
                else if(number_of_subtask==4) task123.new_M_and_F2_FromSpecifiedF1_and_l_OnePart();
                else if(number_of_subtask==5) task123.new_M_and_F1_FromSpecifiedF2_and_l_OnePart();

                if(number_of_subtask!=0 && task123.l!=0.0)
                {
                    qDebug()<<"task123.isTheOptimalPoint_OnePart(task123.l) = "<<task123.isTheOptimalPoint_OnePart(task123.l);
                    info+="\nНовые данные, являющиеся оптимальными:\nl = "+QString("%1").arg(task123.l, 0, 'f', 16)+"\nF_1 = "+QString("%1").arg(task123.F_1, 0, 'f', 16)+"\nF_2 = "+QString("%1").arg(task123.F_2, 0, 'f', 16)+"\nM = "+QString("%1").arg(task123.M, 0, 'f', 16);
                    task123.fillingInVectors_OnePart(x_i, F_i);
                    plotGraphByPoints(x_i, F_i, false);
                }
                else
                {
                    info+="\nДля данной задачи не нашлось решений";
                    qDebug()<<"Для данной задачи не нашлось решений!!!!";
                    QMessageBox::information(this, "Информация о решении", "Для данной задачи не нашлось решений");
                }
            }
        }
        else if(ui->comboBox_TypeOfTask->currentIndex()==1)
        {
            QVector<double> tmp_vector123=task123.decision_TwoParts1();
            info+="\nДля стержня с двумя сечениями(переменное - постоянное) нашлось точек: "+QString::number(tmp_vector123.length())+"\nТочка = "+QString::number(tmp_vector123.value(0));


            if(tmp_vector123.length())
            {
                task123.fillingInVectors_TwoParts1(x_i, F_i, tmp_vector123.value(0));
                plotGraphByPoints(x_i, F_i, false);
                // info+="\nДля двух частей\nx_opt = "+QString::number(task123.x_opt1);
            }
            else
            {
                info+="Данная постановка не удовлетворяет необходимым условиям оптимальности";
                qDebug()<<"Данная постановка не удовлетворяет необходимым условиям оптимальности!!!!";

                QMessageBox::StandardButton reply;
                reply = QMessageBox::question(this, "Test", "Данная постановка не удовлетворяет необходимым условиям оптимальности. Найти оптимальную длину?", QMessageBox::Yes|QMessageBox::No);
                if (reply == QMessageBox::Yes)
                {
                    double tmp12345=task123.searchFor_x_TwoParts1();
                    task123.findOptimalLength_TwoParts1(tmp12345);
                    tmp_vector123=task123.decision_TwoParts1();
                    task123.fillingInVectors_TwoParts1(x_i, F_i, tmp12345);
                    plotGraphByPoints(x_i, F_i, false);
                    info+="\n\nДля двух частей не получилось построить график. После изменений\nl ="+QString("%1").arg(task123.l, 0, 'f', 4)+"\nx_opt = "+QString::number(tmp12345);
                }

            }
        }
        else if(ui->comboBox_TypeOfTask->currentIndex()==2)
        {
            QVector<double> tmp_vector123=task123.decision_TwoParts2();
            info+="\nНашлось точек: "+QString::number(tmp_vector123.length())+"\nТочка = "+QString::number(tmp_vector123.value(0));


            if(tmp_vector123.length())
            {
                task123.fillingInVectors_TwoParts2(x_i, F_i, tmp_vector123.value(0));
                plotGraphByPoints(x_i, F_i, false);
                // info+="\nДля двух частей\nx_opt = "+QString::number(task123.x_opt1);
            }
            else
            {
                qDebug()<<"Данная постановка не удовлетворяет необходимым условиям оптимальности!!!!";
                QMessageBox::StandardButton reply;
                reply = QMessageBox::question(this, "Test", "Данная постановка не удовлетворяет необходимым условиям оптимальности. Найти оптимальную длину?", QMessageBox::Yes|QMessageBox::No);
                if (reply == QMessageBox::Yes)
                {
                    // double x_1_tmp=0.0;
                    task123.findOptimalLength_TwoParts2();
                    task123.fillingInVectors_TwoParts2(x_i, F_i, task123.x_opt1);
                    plotGraphByPoints(x_i, F_i, false);
                    info+="\n\nДля двух частей не получилось построить график. После изменений\nl ="+QString("%1").arg(task123.l, 0, 'f', 3)+"\nx_opt = "+QString("%1").arg(task123.x_opt1, 0, 'f', 3);
                }
            }
        }
        else if(ui->comboBox_TypeOfTask->currentIndex()==3)
        {
            task123.decisionThreeParts();
            info+="\nКол-во найденных решений для трех частей: "+QString::number(task123.x_1_arr.length());
            if(task123.is_task_success==1) QMessageBox::information(this, "Информация о решении", "Для данной задачи не нашлось решений");
            else
            {
                int length_opt_vect=task123.x_1_arr.length();
                for(int i=0; i<length_opt_vect; i++)
                {
                    // info+="\n("+QString::number(task123.x_1_arr.value(i))+", "+QString::number(task123.x_2_arr.value(i))+")";
                    double alpha=task123.findingAlpha(task123.x_1_arr.value(i));
                    info+=QString("\nbeta = %1\nalpha = %2\nx_1 = %3\nx_2 = %4").arg(task123.beta, 0, 'f', 3).arg(alpha, 0, 'f', 3).arg(task123.x_1_arr.value(i), 0, 'f', 3).arg(task123.x_2_arr.value(i), 0, 'f', 3);
                    task123.x_opt1=task123.x_1_arr.value(i);
                    task123.x_opt2=task123.x_2_arr.value(i);
                    x_i.clear();
                    F_i.clear();
                    task123.fillingInVectors3(x_i, F_i, {task123.x_1_arr.value(i), task123.x_2_arr.value(i)});
                    if(i==1) color_graph=Qt::blue;
                    if(i==2) color_graph=Qt::green;

                    plotGraphByPoints(x_i, F_i, false);
                    length_opt_vect=1;  // !!!!!
                }
            }

            task123.checkingTheResultsWithZeroMassAndTheTaskFromTheBook(task123.x_1_arr.value(0), task123.x_2_arr.value(0), task123.findingAlpha(task123.x_1_arr.value(0)));
        }        
        ui->textBrowser_Out->setText(info);
    }
}


void MainWindow::on_comboBox_Material_activated(int index)
{
    ui->comboBox_Material->setToolTip(ui->comboBox_Material->currentText());
    if(index == 0)
    {
        ui->label_MaterialName->setVisible(true);
        ui->label_MaterialDensity->setVisible(true);
        ui->label_MaterialYModule->setVisible(true);

        ui->lineEdit_MaterialName->setVisible(true);
        ui->lineEdit_MaterialDensity->setVisible(true);
        ui->lineEdit_MaterialYModule->setVisible(true);

        ui->pushButton_SaveMaterial->setVisible(true);
    }
    else
    {
        ui->label_MaterialName->setVisible(false);
        ui->label_MaterialDensity->setVisible(false);
        ui->label_MaterialYModule->setVisible(false);

        ui->lineEdit_MaterialName->setVisible(false);
        ui->lineEdit_MaterialDensity->setVisible(false);
        ui->lineEdit_MaterialYModule->setVisible(false);

        ui->pushButton_SaveMaterial->setVisible(false);

        rho=mat123.material_list.value(index-1)->density;
        moduleY=mat123.material_list.value(index-1)->moduleY;
    }
}

void MainWindow::on_pushButton_ColorChange_clicked()
{
    color_graph = QColorDialog::getColor(color_graph, this);
    QString s("background: #"
              + QString(color_graph.red() < 16? "0" : "") + QString::number(color_graph.red(),16)
              + QString(color_graph.green() < 16? "0" : "") + QString::number(color_graph.green(),16)
              + QString(color_graph.blue() < 16? "0" : "") + QString::number(color_graph.blue(),16)
              + ";\nborder-radius: 5px;"
              + "\nborder: 3px solid white;");
    ui->pushButton_ColorChange->setStyleSheet(s);
    ui->pushButton_ColorChange->update();
}
void MainWindow::on_pushButton_Move_clicked(bool checked)
{
    if(checked==true)
    {
        QString s1="background: white;\nborder-radius: 5px;\nborder: 3px solid white;";
        ui->pushButton_Move->setStyleSheet(s1);
        ui->pushButton_Move->update();
        ui->widget_Graph->setInteraction(QCP::iRangeDrag, true);
    }
    else
    {
        QString s1="background: #eaeaea;\nborder-radius: 5px;\nborder: 3px solid #eaeaea;";
        ui->pushButton_Move->setStyleSheet(s1);
        ui->pushButton_Move->update();
        ui->widget_Graph->setInteraction(QCP::iRangeDrag, false);
    }
}


void MainWindow::on_pushButton_ClearGraph_pressed()
{
    QString s1="background: white;\nborder-radius: 5px;\nborder: 3px solid white;";
    ui->pushButton_ClearGraph->setStyleSheet(s1);
    ui->pushButton_ClearGraph->update();
}


void MainWindow::on_pushButton_SaveGraph_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Information"), "E:/Programming/tests/test_save", tr("Text Files (*.txt)"));
    QFile fileOut(fileName);

    if(fileOut.open(QIODevice::WriteOnly | QIODevice::Text))
    { // Если файл успешно открыт для записи в текстовом режиме
        QTextStream writeStream(&fileOut); // Создаем объект класса QTextStream
        // и передаем ему адрес объекта fileOut
        writeStream << printingInformationAboutGraphs(); // Посылаем строку в поток для записи
        fileOut.close(); // Закрываем файл
    }
    else
    {
        qDebug()<<"Ошибка с открытием или записью файла";
    }
}


void MainWindow::on_pushButton_SaveImage_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Image"), "E:/Programming/tests/test_save", tr("Image Files (*.png)"));

    QFile file(fileName);

    if (!file.open(QIODevice::WriteOnly))
    {
        qDebug() << file.errorString();
    } else {
        ui->widget_Graph->savePng(fileName);
    }
}


void MainWindow::on_pushButton_OpenFile_clicked()
{
    OpenTxtFile();
}

void MainWindow::on_pushButton_AddDot_clicked()
{
    double x=ui->doubleSpinBox_InputDot1->value();
    double F=ui->doubleSpinBox_InputDot2->value();
    QVector<double> x_i, F_i;
    x_i.push_back(x);
    F_i.push_back(F);
    // ui->widget_Graph->addGraph();
    // int number=ui->widget_Graph->graphCount();
    ui->widget_Graph->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
    ui->widget_Graph->graph(0)->setLineStyle(QCPGraph::lsNone);

    // ui->widget_Graph->graph(0)->setData(x_i,F_i);
    ui->widget_Graph->graph(0)->addData(x_i,F_i);
    ui->widget_Graph->replot();
    ui->widget_Graph->update();
    qDebug()<<"add point = "<<x;
}


void MainWindow::on_pushButton_Zoom_clicked(bool checked)
{
    if(checked==true)
    {
        QString s1="background: white;\nborder-radius: 5px;\nborder: 3px solid white;";
        ui->pushButton_Zoom->setStyleSheet(s1);
        ui->pushButton_Zoom->update();
        ui->widget_Graph->setInteraction(QCP::iRangeZoom, true);
    }
    else
    {
        QString s1="background: #eaeaea;\nborder-radius: 5px;\nborder: 3px solid #eaeaea;";
        ui->pushButton_Zoom->setStyleSheet(s1);
        ui->pushButton_Zoom->update();
        ui->widget_Graph->setInteraction(QCP::iRangeZoom, false);
    }
}

void MainWindow::on_pushButton_ExampleTask_clicked()
{
    if(ui->comboBox_TypeOfTask->currentIndex()==1)
    {
        ui->lineEdit_l->setText("1.616");
        ui->lineEdit_F_2->setText("1.5");
        ui->lineEdit_F_1->setText("1.0");
        ui->lineEdit_M->setText("100.0");
        ui->lineEdit_omega->setText("5299.9");
    }
    else if(ui->comboBox_TypeOfTask->currentIndex()==2)
    {
        ui->lineEdit_l->setText("3.895");
        ui->lineEdit_F_2->setText("2.0");
        ui->lineEdit_F_1->setText("1.0");
        ui->lineEdit_M->setText("10000.0");
        ui->lineEdit_omega->setText("5299.9");
    }
    else if(ui->comboBox_TypeOfTask->currentIndex()==3)
    {
        ui->lineEdit_l->setText("1.8");
        ui->lineEdit_F_2->setText("2.0");
        ui->lineEdit_F_1->setText("1.0");
        ui->lineEdit_M->setText("100.0");
        ui->lineEdit_omega->setText("5299.9");
    }

    ui->comboBox_Material->setCurrentIndex(1);
    // ui->comboBox_TypeOfTask->setCurrentIndex(1);
    on_comboBox_Material_activated(ui->comboBox_Material->currentIndex());
}


void MainWindow::on_pushButton_SaveMaterial_clicked()
{
    QString material_name_text=ui->lineEdit_MaterialName->text();
    QString material_density_text=ui->lineEdit_MaterialDensity->text();
    QString material_moduleY_text=ui->lineEdit_MaterialYModule->text();

    QString str_if_not_correct="";

    double material_density = stringToNumber(material_density_text, str_if_not_correct, "длины");
    double material_moduleY = stringToNumber(material_moduleY_text, str_if_not_correct, "F_1");


    qDebug()<<"str_if_not_correct = "<<str_if_not_correct;

    if(material_name_text.isEmpty())
    {
        if(str_if_not_correct=="") str_if_not_correct+="названия материала";
        else str_if_not_correct+=", названия материала";
    }
    if(material_density<=0)
    {
        if(str_if_not_correct=="") str_if_not_correct+="плотности материала";
        else str_if_not_correct+=", плотности материала";
    }
    if(material_moduleY<=0)
    {
        if(str_if_not_correct=="") str_if_not_correct+="модуля Юнга материала";
        else str_if_not_correct+=", модуля Юнга материала";
    }


    if(str_if_not_correct!="")
    {
        str_if_not_correct+=". Строки не должны содержать буквы и отрицательные значения.";
        QString info_str="Некорректный ввод данных задачи!\nПроверьте правильность ввода "+str_if_not_correct;
        QMessageBox::warning(this, "Ввод данных", info_str);
    }
    else
    {
        qDebug()<<"Запись в файл материала1";
        qDebug()<<"material_name_text = "<<material_name_text<<"\nmaterial_density = "<<material_density<<"\nmaterial_moduleY = "<<material_moduleY;
        mat123.AddMaterialToFile(material_name_text, material_density, material_moduleY);

        ui->comboBox_Material->addItem(mat123.material_list.value(mat123.material_list.length()-1)->aboutMaterial());
    }
}

void MainWindow::slotForm(int number)
{
    qDebug()<<"\n\n!!!num,ber = "<<number;
    number_of_subtask=number;
    // return number;
}

