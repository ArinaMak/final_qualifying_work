#include "dialogchangeparam.h"
#include "ui_dialogchangeparam.h"

DialogChangeParam::DialogChangeParam(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::DialogChangeParam)
{
    ui->setupUi(this);
}

DialogChangeParam::~DialogChangeParam()
{
    delete ui;
}

void DialogChangeParam::on_buttonBox_accepted()
{
    emit signalForm(ui->comboBox->currentIndex());
}
