#ifndef DIALOGCHANGEPARAM_H
#define DIALOGCHANGEPARAM_H

#include <QDialog>

namespace Ui {
class DialogChangeParam;
}

class DialogChangeParam : public QDialog
{
    Q_OBJECT

public:
    explicit DialogChangeParam(QWidget *parent = nullptr);
    ~DialogChangeParam();

private:
    Ui::DialogChangeParam *ui;

signals:
    void signalForm(int number);

private slots:
    void on_buttonBox_accepted();
};

#endif // DIALOGCHANGEPARAM_H
