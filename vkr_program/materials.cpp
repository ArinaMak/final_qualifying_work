#include "materials.h"
#include <QDir>

Materials::Materials()
{
    OpenMaterialsFile();
}

Materials::~Materials()
{
    int length = material_list.length();
    for(int i=0; i<length; i++)
    {
        delete material_list.value(length - i-1);
        material_list.removeLast();
    }
}

void Materials::OpenMaterialsFile()
{
    // QString fileName = "E:/Programming/fqw/test1_11_05_24/materials.txt";
    QString fileName = ":/rec/text_files/materials2.txt";

    QFile fileIn(fileName);

    if(fileIn.open(QIODevice::ReadOnly))
    {
        QTextStream in(&fileIn);
        in.setCodec("UTF-8");
        int number_of_curr_line=1;

        QString name_material;
        double density=0.0;
        double moduleY=0.0;

        while (!in.atEnd())
        {
            QString line = in.readLine();

            if(number_of_curr_line==1)
            {
                name_material=line;
                qDebug()<<"name_material="<<name_material;
            }
            else if(number_of_curr_line==2)
            {
                density=line.toDouble();
            }
            else if(number_of_curr_line==3)
            {
                moduleY=line.toDouble();
                Material* tmp_mat=new Material(name_material, density, moduleY);
                material_list.push_back(tmp_mat);
                number_of_curr_line=0;
                int tmp_max_length=tmp_mat->aboutMaterial().length();
                if(tmp_max_length>max_length) max_length=tmp_max_length;
            }
            number_of_curr_line++;

        }
        fileIn.close();
    }

    QFile fileIn2("materials123.txt");

    if(fileIn2.open(QIODevice::ReadOnly))
    {
        QTextStream in(&fileIn2);
        in.setCodec("UTF-8");
        int number_of_curr_line=1;

        QString name_material;
        double density=0.0;
        double moduleY=0.0;

        while (!in.atEnd())
        {
            QString line = in.readLine();

            if(number_of_curr_line==1)
            {
                name_material=line;
                qDebug()<<"name_material="<<name_material;
            }
            else if(number_of_curr_line==2)
            {
                density=line.toDouble();
            }
            else if(number_of_curr_line==3)
            {
                moduleY=line.toDouble();
                Material* tmp_mat=new Material(name_material, density, moduleY);
                material_list.push_back(tmp_mat);
                number_of_curr_line=0;
                int tmp_max_length=tmp_mat->aboutMaterial().length();
                if(tmp_max_length>max_length) max_length=tmp_max_length;
            }
            number_of_curr_line++;

        }
        fileIn2.close();
    }

}

void Materials::AddMaterialToFile(QString str, double density, double moduleY)
{
    Material* tmp_mat=new Material(str, density, moduleY);
    material_list.push_back(tmp_mat);

    int tmp_max_length=tmp_mat->aboutMaterial().length();
    if(tmp_max_length>max_length) max_length=tmp_max_length;

    QFile fileOut("materials123.txt");
    if(fileOut.open(QIODevice::WriteOnly | QIODevice::Append))
    { //Если первый файл открыт для чтения, а второй для записи успешн
        QTextStream out(&fileOut); // Создаем объект класса QTextStream
        out.setCodec("UTF-8");
        // и передаем ему адрес объекта fileOut
        out << "\n"+str+"\n"+QString::number(density)+"\n"+QString::number(moduleY); // Посылаем строку в поток для записи

        fileOut.close(); // Закрываем fileout.txt
    }

}
