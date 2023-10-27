/** @file tree.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: A. Bressan.
*/


#pragma once

#include <vector>
#include <iostream> // for debugging and serialization
#include <limits.h>
#include <gsCore/gsDebug.h>


namespace gismo {


// TODO streamline implementation a bit
// add guards values to the iterators:
// values they cannot go through so that
// they can be used to iterate on sub-trees

template <typename NodeId>
struct NodeIdTraits
{
    enum { nullvalue = NodeId(~NodeId(0)) };
};
template <typename NodeId>
struct NodeIdTraits<NodeId*>
{
    enum { nullvalue = NULL };
};

/** A node of a binary tree (see BinaryTree).*/
template <typename Data, typename NodeId>
class TreeNode
{
public:
    static const NodeId no_node=NodeIdTraits<NodeId>::nullvalue;
public:
    NodeId parent;
    NodeId first;
    NodeId second;

    Data   data;
public:
    TreeNode()
        : parent(no_node), first(no_node), second(no_node)
    {}

    bool isFork() const
    {
        return first != no_node || second != no_node;
    }
    bool isLeaf() const
    {
        return first == no_node && second == no_node;
    }
    std::ostream& print(std::ostream &out) const
    {
        out<<"     parent "<<parent<<"\n"
           <<"      first "<<first<<"\n"
           <<"     second "<<second<<"\n"
           <<"       data "<<data;
        (out<<"\n").flush();
        return out;
    }

    operator Data& () {return data;}
    operator const Data& () const {return data;}
};

template <typename Data, typename NodeId>
std::ostream& operator<< (std::ostream& out, const TreeNode<Data,NodeId> &node)
{
    node.print(out);
    return out;
}





template <typename Data>
struct BinaryTree;

template <typename Data>
std::ostream& operator<< (std::ostream& out, const BinaryTree<Data> & tree);


/** Binary tree. The gsDomainMap is its child.*/
template <typename Data>
struct BinaryTree
{
public:
    typedef int NodeId;
    typedef TreeNode<Data,NodeId> Node;
    static const NodeId no_node=NodeIdTraits<NodeId>::nullvalue;
protected:
    std::vector<Node> m_node;
    NodeId            m_freeList;
    NodeId            m_root;
protected:
    NodeId allocateNode()
    {
        if (m_freeList!=no_node)
        {
            NodeId newNode=m_freeList;
            m_freeList=m_node[newNode].parent;
            new (& m_node[newNode]) Node();
            return newNode;
        }
        else
            m_node.resize(m_node.size()+1);
        return m_node.size()-1;
    }


    void freeNode(NodeId node)
    {
        if (m_node[node].parent!=no_node) // Not a root.
        {
            if (node == m_node[m_node[node].parent].first) // We're the first child.
                m_node[m_node[node].parent].first=no_node;
            if (node == m_node[m_node[node].parent].second) // We're the second child.
                m_node[m_node[node].parent].second=no_node;
            m_node[node].parent=m_freeList;
        }
        m_freeList=node;
        m_node[node].data.~Data();

        if (m_node[node].first !=no_node && m_node[node].first != node)
            freeNode(m_node[node].first);
        if (m_node[node].second !=no_node && m_node[node].second != node)
            freeNode(m_node[node].second);
    }
public:
    size_t size () {return m_node.size();}

    size_t depth (const NodeId node)
    {
        size_t result=0;
        NodeId curNode=node;
        for (; m_node[curNode].parent!=no_node; ++result)
            curNode=m_node[curNode].parent;
        return result;
    }

    NodeId append (const NodeId parent, bool first)
    {
        NodeId newNode=allocateNode();
        NodeId *toSet;
        if ( first )
            toSet= &m_node[parent].first;
        else
            toSet= &m_node[parent].second;
        if (*toSet != no_node)
            freeNode(*toSet);
        *toSet=newNode;
        m_node[newNode].parent=parent;
        return newNode;
    }
    /** Cuts away the child specified by firstIsTheOneToCutAway and
     * puts the other on place of node. **/
    void cutBranchAway ( const NodeId node, const bool firstIsTheOneToCutAway )
    {
        // Cutting away the first child.
        if( firstIsTheOneToCutAway && m_node[node].first != no_node)
        {
            // Make the current node ignored.
            if( m_node[node].parent != no_node ) // We're not in the root.
            {
                if( node == m_node[m_node[node].parent].first )
                    m_node[m_node[node].parent].first = m_node[node].second;
                else
                    m_node[m_node[node].parent].second = m_node[node].second;

                m_node[m_node[node].second].parent = m_node[node].parent;
            }
            else // We are in the root: the remaining child becomes the new root.
                setRoot( m_node[node].second );

            // Remove the now isolated branch.
            m_node[node].second = no_node;
            freeNode( node );
        }
        // Cutting away the second child.
        else if ( (!firstIsTheOneToCutAway) && m_node[node].second != no_node )
        {
            if( m_node[node].parent != no_node )
            {
                if( node == m_node[m_node[node].parent].first )
                    m_node[m_node[node].parent].first = m_node[node].first;
                else
                    m_node[m_node[node].parent].second = m_node[node].first;

                m_node[m_node[node].first].parent = m_node[node].parent;
            }
            else
                setRoot( m_node[node].first );

            m_node[node].first = no_node;
            freeNode( node );
        }
        else
        {
            GISMO_ERROR( "You are trying to cut away a non-existing branch!" );
        }
    }
    void remove (const NodeId node)
    {
        freeNode(node);
    }
    void setRoot(const NodeId node)
    {
        if (node==m_root)
            return;
        NodeId parent=m_node[node].parent;
        if (parent != no_node && m_node[parent].first==node)
            m_node[parent].first=no_node;
        if (parent != no_node && m_node[parent].second==node)
            m_node[parent].second=no_node;
        freeNode(m_root);
        m_root=node;
        m_node[node].parent = no_node;
        return;
    }
    void swap(const NodeId node1, const NodeId node2)
    {
        if (node1==no_node)
        {   freeNode(node2); return; }
        if (node2==no_node)
        {   freeNode(node1); return; }
        if (node1==m_root)
        {   setRoot(node2); return; }
        if (node2==m_root)
        {   setRoot(node1); return; }

        NodeId par1=m_node[node1].parent;
        NodeId par2=m_node[node2].parent;

        m_node[node2].parent=par1;
        m_node[node1].parent=par2;

        if ( (par1!=no_node) && (m_node[par1].first==node1) )
            m_node[par1].first=node2;
        if ( (par1!=no_node) && (m_node[par1].second==node1) )
            m_node[par1].second=node2;
        if ( (par2!=no_node) && (m_node[par2].first==node2) )
            m_node[par2].first=node1;
        if ( (par2!=no_node) && (m_node[par2].second==node2) )
            m_node[par2].second=node1;
    }
    void detach (const NodeId node)
    {
        if (node==no_node)
        {
            gsWarn<<"HELLO!!!\n";
            return;
        }
        NodeId parent=m_node[node].parent;
        if( parent != no_node && m_node[parent].first==node)
            m_node[parent].first=no_node;
        if( parent != no_node && m_node[parent].second==node)
            m_node[parent].second=no_node;
        m_node[node].parent=no_node;
        if (m_root==node)
            m_root=no_node;
    }

    NodeId root() const
    {
        return m_root;
    }
    
    Node& operator[] (const NodeId node)
    {
        return m_node[node];
    }

    const Node& operator[] (const NodeId node) const
    {
        return m_node[node];
    }
protected:
    struct iterate_fork_last
    {
        static inline NodeId first(const BinaryTree &tree, NodeId root)
        {
            while (true)
            {
                while (tree.m_node[root].first != no_node )
                    root=tree.m_node[root].first;
                if ( tree.m_node[root].second != no_node )
                    root=tree.m_node[root].second;
                else
                    return root;
            }
        }
        static inline NodeId next (const BinaryTree &tree, NodeId cur)
        {
            NodeId parent=tree.m_node[cur].parent;
            if ( parent == no_node || tree.m_node[parent].second==cur || tree.m_node[parent].second==no_node )
                return parent;

            cur=tree.m_node[parent].second;
            while (true)
            {
                while (tree.m_node[cur].first != no_node)
                    cur=tree.m_node[cur].first;
                if( tree.m_node[cur].second != no_node)
                    cur=tree.m_node[cur].second;
                else return cur;
            }
            return tree.m_node[parent].second;
        }
        static inline NodeId prev (const BinaryTree &tree, NodeId cur)
        {
            if ( tree.m_node[cur].second != no_node )
                return tree.m_node[cur].second;
            else if ( tree.m_node[cur].first != no_node )
                return tree.m_node[cur].first;
            else while (true)
            {
                NodeId parent=tree.m_node[cur].parent;
                if ( parent == no_node )
                    return no_node;
                else if (tree.m_node[parent].first != no_node && tree.m_node[parent].first != cur)
                    return tree.m_node[parent].first;
                else
                    cur=parent;
            }
        }
        static inline NodeId last (const BinaryTree &tree, NodeId root)
        {
            return root;
        }
    };
    struct iterate_fork_first
    {
        static inline NodeId last(const BinaryTree &tree, NodeId root)
        {
            while (true)
            {
                while (tree.m_node[root].second != no_node )
                    root=tree.m_node[root].second;
                if ( tree.m_node[root].first != no_node )
                    root=tree.m_node[root].first;
                else
                    return root;
            }
        }
        static inline NodeId prev (const BinaryTree &tree, NodeId cur)
        {
            NodeId parent=tree.m_node[cur].parent;
            if ( parent == no_node || cur == tree.m_node[parent].first || tree.m_node[parent].first==no_node)
                return parent;

            cur=tree.m_node[parent].first;
            while (true)
            {
                while (tree.m_node[cur].second != no_node)
                    cur=tree.m_node[cur].second;
                if( tree.m_node[cur].first != no_node)
                    cur=tree.m_node[cur].first;
                else return cur;
            }
            return tree.m_node[parent].first;
        }
        static inline NodeId next (const BinaryTree &tree, NodeId cur)
        {
            if ( tree.m_node[cur].first != no_node )
                return tree.m_node[cur].first;
            else if ( tree.m_node[cur].second != no_node )
                return tree.m_node[cur].second;
            else while (true)
            {
                NodeId parent=tree.m_node[cur].parent;
                if ( parent == no_node )
                    return no_node;
                else if (tree.m_node[parent].second != no_node && tree.m_node[parent].second != cur)
                    return tree.m_node[parent].second;
                else
                    cur=parent;
            }
        }
        static inline NodeId first (const BinaryTree & /*tree*/,
                                          NodeId       root)
        {
            return root;
        }
    };
    struct iterate_fork_middle
    {
        static inline NodeId first(const BinaryTree &tree, NodeId root)
        {
            while (true)
            {
                while (tree.m_node[root].first != no_node )
                    root=tree.m_node[root].first;
                return root;
            }
        }
        static inline NodeId next (const BinaryTree &tree, NodeId cur)
        {
            if ( tree.m_node[cur].second != no_node)
            {
                cur=tree.m_node[cur].second;
                while ( tree.m_node[cur].first != no_node)
                    cur=tree.m_node[cur].first;
                return cur;
            }
            else while (true)
            {
                NodeId parent=tree.m_node[cur].parent;
                if ( parent==no_node || tree.m_node[parent].first==cur)
                    return parent;
                else
                    cur=parent;
            }
        }
        static inline NodeId prev (const BinaryTree &tree, NodeId cur)
        {
            if ( tree.m_node[cur].first != no_node)
            {
                cur=tree.m_node[cur].first;
                while ( tree.m_node[cur].second != no_node)
                    cur=tree.m_node[cur].second;
                return cur;
            }
            else while (true)
            {
                NodeId parent=tree.m_node[cur].parent;
                if ( parent==no_node || tree.m_node[parent].second==cur)
                    return parent;
                else
                    cur=parent;
            }
        }
        static inline NodeId last(const BinaryTree &tree, NodeId root)
        {
            while (true)
            {
                while (tree.m_node[root].second != no_node )
                    root=tree.m_node[root].second;
                return root;
            }
        }
    };
    template <typename base>
    struct revert
    {
        static inline NodeId first(const BinaryTree &tree, NodeId root)     {return base::last(tree, root);}
        static inline NodeId last(const BinaryTree &tree, NodeId root)      {return base::first(tree, root);}
        static inline NodeId next (const BinaryTree &tree, NodeId cur)      {return base::prev(tree, cur);}
        static inline NodeId prev (const BinaryTree &tree, NodeId cur)      {return base::next(tree, cur);}
    };
    template <typename base>
    struct leaf
    {
        static inline NodeId first(const BinaryTree &tree, NodeId root)
        {
            NodeId tmp=base::first(tree, root);
            while (tree[tmp].isFork())
                tmp=base::next(tree,tmp);
            return tmp;
        }
        static inline NodeId last(const BinaryTree &tree, NodeId root)
        {
            NodeId tmp=base::last(tree, root);
            while (tree[tmp].isFork())
                tmp=base::prev(tree,tmp);
            return tmp;
        }
        static inline NodeId next (const BinaryTree &tree, NodeId cur)
        {
            cur=base::next(tree,cur);
            while (cur!=no_node && tree[cur].isFork())
                cur=base::next(tree,cur);
            return cur;
        }
        static inline NodeId prev (const BinaryTree &tree, NodeId cur)
        {
            cur=base::prev(tree,cur);
            while (cur!=no_node && tree[cur].isFork())
                cur=base::prev(tree,cur);
            return cur;
        }
    };


    template <typename order>
    class const_iterator_base;
    template <typename base>
    class iterator_base {
    public:
        // template <typename otherBase>
        // friend class iterator_base;
    protected:
        BinaryTree    *m_tree;
        NodeId         m_node;
    public:
        typedef base order;

        typedef Node  value_type;
        typedef Node& reference;
        typedef Node* pointer;
        typedef std::bidirectional_iterator_tag iterator_category; //or another tag

        iterator_base() : m_tree(NULL), m_node(no_node) {}
        iterator_base(BinaryTree &tree, NodeId node=no_node) : m_tree(&tree), m_node(node) {}

        // template <typename otherBase>
        // iterator_base(const iterator_base<otherBase>& other) :
        // m_tree(other.m_tree),  m_node(other.m_node){}

        ~iterator_base() {}

        operator pointer() {return &(m_tree->m_node[m_node]);}
        operator NodeId () {return m_node;}

        iterator_base& operator=(const iterator_base& other) {m_node=other.m_node;m_tree=other.m_tree; return *this;}
        bool operator==(const iterator_base& other) const     {return  m_node==other.m_node&&m_tree==other.m_tree;}
        bool operator!=(const iterator_base& other) const     {return  m_node!=other.m_node||m_tree!=other.m_tree;}

        iterator_base& operator++()                     {m_node=base::next(*m_tree,m_node); return *this;}
        iterator_base& operator--()                     {m_node=base::prev(*m_tree,m_node); return *this;}
        iterator_base& operator+=(size_t a)             {for (int i=0; i<a;++i) m_node=base::next(*m_tree,m_node); return *this;}
        iterator_base  operator+ (size_t a)   const     {iterator_base other(*this); return other+=a;}
        iterator_base& operator-=(size_t a)             {for (int i=0; i<a;++i) m_node=base::prev(*m_tree,m_node); return *this;}
        iterator_base  operator- (size_t a)   const     {iterator_base other(*this); return other-=a;}


        reference operator*()  const               {return m_tree->m_node[m_node];}
        pointer   operator->() const               {return &(m_tree->m_node[m_node]);}

        typedef iterator_base<revert<order> > reverse_iterator;

        template <typename other>
        operator iterator_base<other>() const {return iterator_base<other>(*m_tree,m_node);}
        template <typename other>
        operator const_iterator_base<other>() const {return const_iterator_base<other>(*m_tree,m_node);}
    };
    template <typename order>
    class const_iterator_base {
    protected:
        const BinaryTree    *m_tree;
        NodeId         m_node;
    public:
        typedef Node        value_type;
        typedef Node const& reference;
        typedef Node const* pointer;
        typedef std::bidirectional_iterator_tag iterator_category; //or another tag

        const_iterator_base() : m_tree(NULL), m_node(no_node) {}
        const_iterator_base(const BinaryTree &tree, NodeId node=no_node) : m_tree(&tree), m_node(node) {}
        const_iterator_base(const iterator_base<order> & other) :m_tree(other.m_tree),  m_node(other.m_node){}
        ~const_iterator_base() {}

        operator pointer() {return &(m_tree->m_node[m_node]);}
        operator NodeId () {return m_node;} // this allows to circumvent const but it is inevitable

        const_iterator_base& operator=(const const_iterator_base& other) {m_node=other.m_node; return *this;}
        bool operator==(const const_iterator_base& other) const     {return  m_node==other.m_node;}
        bool operator!=(const const_iterator_base& other) const     {return  m_node!=other.m_node;}

        const_iterator_base& operator++()                     {m_node=order::next(*m_tree,m_node); return *this;}
        const_iterator_base& operator--()                     {m_node=order::prev(*m_tree,m_node); return *this;}
        const_iterator_base& operator+=(size_t a)             {for (int i=0; i<a;++i) m_node=order::next(*m_tree,m_node); return *this;}
        const_iterator_base  operator+ (size_t a)   const     {const_iterator_base other(*this); return other+=a;}
        const_iterator_base& operator-=(size_t a)             {for (int i=0; i<a;++i) m_node=order::prev(*m_tree,m_node); return *this;}
        const_iterator_base  operator- (size_t a)   const     {const_iterator_base other(*this); return other-=a;}

        reference operator*()  const               {return m_tree->m_node[m_node];}
        pointer   operator->() const               {return &(m_tree->m_node[m_node]);}

        typedef const_iterator_base<revert<order> > reverse_iterator;

        template <typename other>
        operator const_iterator_base<other>() const {return const_iterator_base<other>(*m_tree,m_node);}
    };
public:
    // orders for nodes
    typedef iterate_fork_first          all_NFS;
    typedef revert<iterate_fork_first>  all_SFN;
    typedef iterate_fork_last           all_FSN;
    typedef revert<iterate_fork_last>   all_NSF;
    typedef iterate_fork_middle         all_FNS;
    typedef revert<iterate_fork_middle> all_SNF;

    // orders for leafs
    typedef leaf<all_NFS>               leaf_FS;
    typedef leaf<all_SFN>               leaf_SF;

    // node iterators
    typedef iterator_base<all_NFS>    node_iterator_NFS;
    typedef iterator_base<all_SFN>    node_iterator_SFN;
    typedef iterator_base<all_FSN>    node_iterator_FSN;
    typedef iterator_base<all_NSF>    node_iterator_NSF;
    typedef iterator_base<all_FNS>    node_iterator_FNS;
    typedef iterator_base<all_SNF>    node_iterator_SNF;


    // leaf iterators
    typedef iterator_base<leaf_FS>   leaf_iterator_FS;
    typedef iterator_base<leaf_SF>   leaf_iterator_SF;

    // const node iterators
    typedef const_iterator_base<all_NFS>    const_node_iterator_NFS;
    typedef const_iterator_base<all_SFN>    const_node_iterator_SFN;
    typedef const_iterator_base<all_FSN>    const_node_iterator_FSN;
    typedef const_iterator_base<all_NSF>    const_node_iterator_NSF;
    typedef const_iterator_base<all_FNS>    const_node_iterator_FNS;
    typedef const_iterator_base<all_SNF>    const_node_iterator_SNF;


    // const leaf iterators
    typedef const_iterator_base<leaf_FS>   const_leaf_iterator_FS;
    typedef const_iterator_base<leaf_SF>   const_leaf_iterator_SF;

    // functions
    NodeId getRoot () const
    {
        return m_root;
    }

    template <typename order>
    iterator_base<order> begin()  { return iterator_base<order>(*this, order::first(*this, m_root)); }
    template <typename order>
    iterator_base<order> last()   { return iterator_base<order>(*this, order::last(*this, m_root)); }
    template <typename order>
    iterator_base<order> end()    { return iterator_base<order>(*this,no_node); }

    template <typename order>
    const_iterator_base<order> begin() const { return const_iterator_base<order>(*this, order::first(*this, m_root)); }
    template <typename order>
    const_iterator_base<order> last()  const { return const_iterator_base<order>(*this, order::last(*this, m_root)); }
    template <typename order>
    const_iterator_base<order> end()   const { return const_iterator_base<order>(*this,no_node); }

public:
    BinaryTree ()
    {
        m_freeList=no_node;
        m_root=0;
        m_node.resize(1);
    }

    // warning : no deep copy
    BinaryTree ( const BinaryTree &other, bool compact)
    {
        if (!compact)
        {
            m_freeList=other.m_freeList;
            m_root=other.m_root;
            m_node=other.m_node;
        }
        else
        {
            const_node_iterator_NFS oldNode = other.template begin<all_NFS>();
            node_iterator_NFS newNode = begin<all_NFS>();
            m_freeList=no_node;
            // preconstruct root
            m_root=0;
            m_node.resize(1);
            // init iterations
            while(oldNode)
            {
                newNode->data=oldNode->data;
                NodeId fst;
                NodeId snd;
                fst = (oldNode->first  != no_node)? append(newNode,true)  : no_node;
                snd = (oldNode->second != no_node)? append(newNode,false) : no_node;
                ++newNode;
                ++oldNode;
            }
        }
    }

    friend std::ostream& operator<< <> (std::ostream& out, const BinaryTree & tree);
};

template <typename Data>
std::ostream& operator<< (std::ostream &out, const BinaryTree<Data> &tree);

}
