#include "material.h"

Material::Material(QString sname, double sdensity, double smoduleY)
{
    name=sname;
    density=sdensity;
    moduleY=smoduleY;
}

QString Material::aboutMaterial()
{
    qDebug()<<"name = "<<name<<"\tdensity = "<<density<<"\tmoduleY = "<<moduleY;
    QString tmp;
    tmp=name+" ( rho="+QString::number(density)+QString::fromUtf8("кг/м\u{00b3}")+", E="+QString::number(moduleY)+")";
    return tmp;
}
