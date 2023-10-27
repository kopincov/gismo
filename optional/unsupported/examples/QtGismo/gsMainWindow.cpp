#include "gsMainWindow.h"

#include "ui_gsMainWindow.h"
#include <gsHSplines/gsQuadTree.h>
using namespace gismo;

gsMainWindow::gsMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::gsMainWindow)
{
    ui->setupUi(this);
}

gsMainWindow::~gsMainWindow()
{
    delete ui;
}

void gsMainWindow::paintEvent(QPaintEvent *e)
{
    quadtree<unsigned int> tree;
    QPoint p1;
    QPoint p2;
    gsVector<unsigned int> i1;
    gsVector<unsigned int> i2;


    //tree.set_knots(0, 21, 20, 3, 1);
    gsKnotVector<> KV (0, 23, 21, 3, 1);
    //tree.knotsX  = KV;
    gsKnotVector<> Kw (0, 16, 15, 3, 1);
    //tree.knotsY = Kw;
   // tree.knotsX(0, 21, 20, 3, 1);
    //tree.knotsY(0, 16, 15, 3, 1);
    //tree.set_unique_knots();
    tree.set_root(64,64);
    vector<real_t> q;
    i1[0] = 5;
    i1[1] = 7;
    i2[0] = 10;
    i2[1] = 12;
    tree.insert_box(i1,i2,tree.root,2);
    //q = tree.draw_tree(tree.root,16,true);
    QPainter painter(this);
    //painter.drawLine(10, 10 , 100, 100);
    //QPen pen(Qt::black);
   // pen.setColor(Qt::green);
   // painter.setPen(pen);
    for(unsigned i = 0; i < q.size();i++){
        p1.setX(q[i]);
        p1.setY(q[i+1]+15);
        p2.setX(q[i+2]);
        p2.setY(q[i+3]+15);
        i = i+3;
        painter.drawLine(p1,p2);
    }
    if (tree.query1(i1,i2,2,tree.root)){
        cout<< "Q1 true"<< endl;
    }else{
        cout<< "Q1 false"<< endl;
    }
    i1[0] = 4;
    if (tree.query2(i1,i2,2,tree.root)){
        cout<< "Q1 true"<< endl;
    }else{
        cout<< "Q1 false"<< endl;
    }
    i1[0] = 6;
    int p = tree.query3(i1,i2,tree.root);
    cout<< "minimal level: "<< p << endl;

}

