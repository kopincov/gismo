#ifndef GSMAINWINDOW_H
#define GSMAINWINDOW_H

#include <QMainWindow>
#include <QtGui>
#include <QtCore>
#include <gsHSplines/gsQuadTree.h>

namespace Ui {
    class gsMainWindow;
}

class gsMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit gsMainWindow(QWidget *parent = 0);
    ~gsMainWindow();

private:
    Ui::gsMainWindow *ui;
    QGraphicsScene *scene;

protected:
    void paintEvent(QPaintEvent *e);

};

#endif // GSMAINWINDOW_H
