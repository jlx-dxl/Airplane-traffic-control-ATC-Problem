/*
 * 	$Id: ControlPanel.java,v 1.2 2007/11/13 02:49:55 johnc Exp $
 * 
 * 	Programming and Problem Solving
 *  Copyright (c) 2007 The Trustees of Columbia University
 */
package airplane.sim.ui;

import java.awt.Dimension;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public final class ControlPanel extends JPanel{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected JButton begin;
	protected JButton step;
	protected JButton play;
	protected JButton pause; 
	protected JButton stop;
	protected JTextField roundText;
	protected JTextField powerText;
	protected JTextField delayText;
	protected JButton boardPanel;
	
	protected ControlPanel(){
		begin = new JButton("Begin New Game");
		begin.setName("Begin");
		begin.setActionCommand("Begin");
		begin.setEnabled(true);
		
		step = new JButton("Step");
		step.setName("Step");
		step.setEnabled(false);
		step.setActionCommand("Step");
		
		play = new JButton("Play");
		play.setName("Play");
		play.setEnabled(false);
		play.setActionCommand("Play");
		
		pause = new JButton("Pause");
		pause.setName("Pause");
		pause.setEnabled(false);
		pause.setActionCommand("Pause");
		
		stop = new JButton("Resign");
		stop.setName("Stop");
		stop.setEnabled(false);
		stop.setActionCommand("Stop");
		
		/*
		tournament = new JButton("Run Tournament");
		tournament.setName("Tournament");
		tournament.setEnabled(true);
		tournament.setActionCommand("Tournament");
		
		boardPanel = new JButton("2nd Display");
		boardPanel.setName("BoardFrame");
		boardPanel.setEnabled(true);
		boardPanel.setActionCommand("BoardFrame");
		 */
		
		add(begin);
		add(step);
		add(play);
		add(pause);
		add(stop);
		//add(tournament);
		//add(boardPanel);
		
		JLabel label = new JLabel("Time: ");
		roundText = new JTextField();
		roundText.setEnabled(false);
		roundText.setEditable(false);
		roundText.setPreferredSize(new Dimension(50, 25));
		add(label);
		add(roundText);
		
		label = new JLabel("Power: ");
		powerText = new JTextField();
		powerText.setEnabled(false);
		powerText.setEditable(false);
		powerText.setPreferredSize(new Dimension(50, 25));
		add(label);
		add(powerText);

		label = new JLabel("Delay: ");
		delayText = new JTextField();
		delayText.setEnabled(false);
		delayText.setEditable(false);
		delayText.setPreferredSize(new Dimension(50, 25));
		add(label);
		add(delayText);
	}
	protected void addListener(ActionListener a){
		begin.addActionListener(a);
		pause.addActionListener(a);
		play.addActionListener(a);
		step.addActionListener(a);
		stop.addActionListener(a);
		//tournament.addActionListener(a);
		//boardPanel.addActionListener(a);
	}
}
