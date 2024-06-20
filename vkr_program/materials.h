#ifndef MATERIALS_H
#define MATERIALS_H

#include <QVector>
#include <material.h>
#include <QFile>
class Materials
{
public:
    QVector<Material*> material_list;
    int max_length=0;

    Materials();
    ~Materials();
    // QString aboutMaterial(int i);
    void OpenMaterialsFile();
    void AddMaterialToFile(QString str, double density, double moduleY);
};

#endif // MATERIALS_H
