#ifndef MATERIAL_H
#define MATERIAL_H

#include <QString>
#include <QDebug>

class Material
{
public:
    QString name;
    double density;
    double moduleY;

    Material(QString sname, double sdensity, double smoduleY);
    QString aboutMaterial();
};

#endif // MATERIAL_H
